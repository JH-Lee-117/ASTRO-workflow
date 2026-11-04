#!/usr/bin/env python3
"""
GFAP Phantom Benchmark: LF-calibrated v3 renderer, enhanced proximal process width/intensity.
Inputs (relative to --project-root):
  data/swc_final_low18_medium12/NMO_<ID>.swc
  data/metadata/selected_30_low18_medium12.csv

Outputs (new versioned folder; originals untouched):
  results/phantom_v3_LF/skeleton                       : 30 TIFF
  results/phantom_v3_LF/gt_binary                      : 30 TIFF
  results/phantom_v3_LF/grayscale_clean                : 30 TIFF
  results/phantom_v3_LF/grayscale_moderate             : 30 TIFF
  results/phantom_v3_LF/qc_30.csv
  results/phantom_v3_LF/preview_clean_moderate_30.png
  results/phantom_v3_LF/rendering_config.json

IMPORTANT:
- Coordinate sampling, image size and seed are read from the explicit --config JSON.
- Changing sampling does not rescale source SWC coordinates. Pixel-based model widths
  retain their pixel units, so changing sampling also changes their physical widths.
- Diameter unavailable or unreliable in some SWCs; modeled widths are not measured radii.
- GT masks derive independently from SWC geometry, NEVER from grayscale thresholding.
- Grayscale-to-binary is to be performed with the user's ORIGINAL Fiji workflow.
- Fiji-based reference and prediction Sholl must use same center and 2:2:40 um radii.
- No Python Sholl values are generated: both arms are measured by Fiji.
- Clean and Moderate share precisely the same underlying morphology and GT mask.
"""


from pathlib import Path
import argparse as _argparse
_parser = _argparse.ArgumentParser(description="Render GFAP phantom images from selected SWCs.")
_parser.add_argument('--project-root', type=Path, required=True, help="Analysis project containing data/ and results/.")
_parser.add_argument("--config", type=Path, required=True, help="Input JSON with pixel_size_um, image_size_px and seed; existing rendering_config.json is supported.")
_args = _parser.parse_args()
_project = _args.project_root.resolve()
from collections import defaultdict
from hashlib import sha256
import json
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.draw import line, disk, ellipse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


import json
import math

def load_settings(path):
    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate configuration key: {key}")
            result[key] = value
        return result
    def invalid_constant(value):
        raise ValueError(f"Non-finite JSON value: {value}")
    settings = json.loads(path.read_text(encoding="utf-8"),
                          object_pairs_hook=unique_object, parse_constant=invalid_constant)
    if not isinstance(settings, dict):
        raise ValueError("Configuration must be a JSON object.")
    for key in ("pixel_size_um", "image_size_px", "seed"):
        if key not in settings:
            raise ValueError(f"Configuration missing {key}")
        if type(settings[key]) not in (int, float) or not math.isfinite(settings[key]):
            raise ValueError(f"Configuration {key} must be a finite number.")
    if settings["pixel_size_um"] <= 0:
        raise ValueError("pixel_size_um must be positive.")
    for key, lower in (("image_size_px", 2), ("seed", 0)):
        if settings[key] != int(settings[key]) or settings[key] < lower:
            raise ValueError(f"{key} must be an integer >= {lower}.")
        settings[key] = int(settings[key])
    return settings

SETTINGS = load_settings(_args.config)

def check_selection_settings(table):
    expected = {"pixel_size_um": SETTINGS["pixel_size_um"],
                "image_size_px": SETTINGS["image_size_px"],
                "selection_seed": SETTINGS["seed"]}
    for key, value in expected.items():
        if key not in table.columns:
            raise RuntimeError(f"Existing selection lacks {key}; cannot verify configuration.")
        values = pd.to_numeric(table[key], errors="raise")
        if values.isna().any() or not np.isfinite(values).all() or not (values == value).all():
            raise RuntimeError(f"Selection {key} differs from --config; use a separate project.")

N = SETTINGS["image_size_px"]
PX = float(SETTINGS["pixel_size_um"])
SEED = SETTINGS["seed"]
CENTER = N // 2

ROOT = _project
SOURCES = ROOT / "data/swc_final_low18_medium12"
MANIFEST = ROOT / "data/metadata/selected_30_low18_medium12.csv"
OUT = ROOT / "results/phantom_v3_LF"

CONFIG = {
    "pixel_size_um": PX, "image_size_px": N, "seed": SEED,
    "soma_ellipse_radii_px": [24.0, 18.0],
    "soma_scale_range": [0.85, 1.15],
    "calibration_note": "LF pooled soma area median 31.853 um2 (n=55); method-defined macro area, not anatomical soma boundary",
    "process_radius_px": "cell_factor*(1.5 + 9.0*exp(-path_length_um/10))",
    "cell_factor_range": [0.9, 1.1],
    "clean": {"psf_sigma_px": 0.8, "gain": 1.0, "background": 18, "bg_variation": 2,
              "photon_factor": 4.0, "read_noise_sd": 2.0},
    "moderate": {"psf_sigma_px": 1.3, "gain": 0.76, "background": 24, "bg_variation": 5,
                 "photon_factor": 1.0, "read_noise_sd": 5.0},
    "validation_sholl": "Fiji reference/test masks only, same soma coordinate, 2-40um step2; independent Fiji GT measurement not performed here",
    "proximal_brightness": "65*exp(-(radial_distance_um/10)^2), with baseline 145 and soma contribution 90, pre-noise",
    "revision_note": "v3: LF soma ellipse same as v2; increased proximal process width and brightness; distal width nearly unchanged",
    "geometry": "Native SWC XY in micrometers, no scaling, soma center (200,200)",
}
# Numerical entries control rendering; descriptive model strings document fixed equations.
unknown = set(SETTINGS) - set(CONFIG) - {"selection_margin_um"}
if unknown:
    raise ValueError(f"Unrecognized configuration keys: {sorted(unknown)}")
for key in ("process_radius_px", "proximal_brightness"):
    if key in SETTINGS and SETTINGS[key] != CONFIG[key]:
        raise ValueError(f"{key} documents a fixed equation; changing its text cannot change the model.")
CONFIG.update({key: value for key, value in SETTINGS.items() if key != "selection_margin_um"})
CONFIG["geometry"] = f"Native SWC XY in micrometers, no scaling, soma center ({CENTER},{CENTER})"

def positive_number(value, name, allow_zero=False):
    if (type(value) not in (int, float) or not math.isfinite(value)
            or (value < 0 if allow_zero else value <= 0)):
        raise ValueError(f"Invalid numerical configuration: {name}")

for key in ("soma_ellipse_radii_px", "soma_scale_range", "cell_factor_range"):
    pair = CONFIG[key]
    if not isinstance(pair, list) or len(pair) != 2:
        raise ValueError(f"{key} requires two numerical values.")
    for value in pair:
        positive_number(value, key)
    if key.endswith("range") and pair[0] > pair[1]:
        raise ValueError(f"Invalid range: {key}")
for condition in ("clean", "moderate"):
    values = CONFIG[condition]
    required = {"psf_sigma_px", "gain", "background", "bg_variation", "photon_factor", "read_noise_sd"}
    if not isinstance(values, dict) or set(values) != required:
        raise ValueError(f"{condition} must contain exactly {sorted(required)}")
    for key, value in values.items():
        positive_number(value, condition + "." + key, allow_zero=(key != "photon_factor"))
SOMA_AXES_PX = tuple(CONFIG["soma_ellipse_radii_px"])

def save_tif(path, data):
    if path.exists():
        with tifffile.TiffFile(path) as tif:
            old = tif.asarray()
            metadata = tif.imagej_metadata or {}
            if (old.shape != data.shape or old.dtype != data.dtype
                    or tif.series[0].axes != "YX"):
                raise RuntimeError(f"Existing image format differs; not overwriting: {path}")
            if metadata.get("unit") not in {"um", "µm", "μm"}:
                raise RuntimeError(f"Existing image unit differs; not overwriting: {path}")
            for name in ("XResolution", "YResolution"):
                tag = tif.pages[0].tags.get(name)
                if tag is None:
                    raise RuntimeError(f"Existing image lacks {name}: {path}")
                numerator, denominator = tag.value
                if denominator == 0 or not np.isclose(
                        numerator / denominator, 1 / PX, rtol=1e-6, atol=0):
                    raise RuntimeError(f"Existing image calibration differs; not overwriting: {path}")
            resolution_unit = tif.pages[0].tags.get("ResolutionUnit")
            if resolution_unit is None or int(resolution_unit.value) != 1:
                raise RuntimeError(f"Existing image resolution unit differs; not overwriting: {path}")
        if not np.array_equal(old, data):
            raise RuntimeError(f"Existing image differs; not overwriting: {path}")
        return
    tifffile.imwrite(path, data, imagej=True, metadata={"unit": "um", "axes": "YX"},
                      resolution=(1 / PX, 1 / PX))

def read_swc(path):
    swc = np.atleast_2d(np.loadtxt(path, comments="#"))
    if swc.shape[1] != 7 or not np.isfinite(swc).all():
        raise ValueError(f"Invalid SWC: {path}")
    ids = swc[:, 0].astype(int)
    parents = swc[:, 6].astype(int)
    lookup = {int(n): i for i, n in enumerate(ids)}
    roots = np.flatnonzero(parents == -1)
    if len(roots) != 1 or len(lookup) != len(ids):
        raise ValueError(f"Invalid IDs or roots: {path}")
    children = defaultdict(list)
    edges = []
    for i, p in enumerate(parents):
        if p != -1:
            if int(p) not in lookup:
                raise ValueError(f"Missing parent: {path}")
            j = lookup[int(p)]
            children[j].append(i)
            edges.append((j, i))
    dist = np.full(len(ids), np.nan)
    dist[roots[0]] = 0.0
    stack = [int(roots[0])]
    visited = set()
    while stack:
        j = stack.pop()
        if j in visited:
            raise ValueError(f"Cycle: {path}")
        visited.add(j)
        for i in children[j]:
            dist[i] = dist[j] + np.linalg.norm(swc[i, 2:5] - swc[j, 2:5])
            stack.append(i)
    if len(visited) != len(ids):
        raise ValueError(f"Disconnected: {path}")
    soma_nodes = swc[swc[:, 1].astype(int) == 1, 2:4]
    if not len(soma_nodes):
        raise ValueError(f"No soma: {path}")
    soma = soma_nodes.mean(axis=0)
    xy = swc[:, 2:4] - soma
    cols = np.rint(CENTER + xy[:, 0] / PX).astype(int)
    rows = np.rint(CENTER - xy[:, 1] / PX).astype(int)
    if (cols.min() < 0 or cols.max() >= N or rows.min() < 0 or rows.max() >= N):
        raise ValueError(f"SWC centerlines outside FOV: {path}")
    return swc, edges, dist, soma, rows, cols, xy

def main():
    selected = pd.read_csv(MANIFEST)
    check_selection_settings(selected)
    if len(selected) != 30 or selected.neuron_id.nunique() != 30:
        raise RuntimeError("Expected 30 unique selected SWC IDs.")
    if "complexity_group" not in selected.columns:
        raise RuntimeError("Manifest missing complexity_group; cannot verify Low/Medium-only selection.")
    group_counts = selected["complexity_group"].value_counts().to_dict()
    if group_counts != {"Low": 18, "Medium": 12}:
        raise RuntimeError(f"Invalid selection; expected Low=18, Medium=12, High=0; got {group_counts}")
    print("SELECTION VERIFIED: Low 18, Medium 12, High 0")
    for row in selected.itertuples(index=False):
        path = SOURCES / f"NMO_{int(row.neuron_id)}.swc"
        if not path.is_file():
            raise FileNotFoundError(path)
        if hasattr(row, "sha256") and pd.notna(row.sha256):
            if sha256(path.read_bytes()).hexdigest() != str(row.sha256):
                raise RuntimeError(f"Source checksum mismatch: {path}")

    conf_path = OUT / "rendering_config.json"
    if conf_path.exists() and json.loads(conf_path.read_text()) != CONFIG:
        raise RuntimeError("Rendering config changed; use a new output folder/version.")
    for name in ("skeleton", "gt_binary", "grayscale_clean", "grayscale_moderate"):
        (OUT / name).mkdir(parents=True, exist_ok=True)
    conf_path.write_text(json.dumps(CONFIG, indent=2, ensure_ascii=False))

    qc, preview = [], []
    for row in selected.itertuples(index=False):
        nid = int(row.neuron_id)
        path = SOURCES / f"NMO_{nid}.swc"
        swc, edges, dist, soma, rr, cc, xy = read_swc(path)
        rng_geometry = np.random.default_rng(SEED + nid)
        width_factor = float(rng_geometry.uniform(*CONFIG["cell_factor_range"]))
        radius_nodes = width_factor*(1.5 + 9.0*np.exp(-dist/10.0))

        center_radius = np.zeros((N, N), dtype=np.float32)
        center_radius[rr, cc] = radius_nodes
        for parent, child in edges:
            ys, xs = line(rr[parent], cc[parent], rr[child], cc[child])
            t = np.linspace(0, 1, len(ys))
            r = radius_nodes[parent]*(1 - t) + radius_nodes[child]*t
            np.maximum.at(center_radius, (ys, xs), r)

        skeleton = (center_radius > 0).astype(np.uint8) * 255
        d, nearest = ndi.distance_transform_edt(center_radius == 0, return_indices=True)
        gt = d <= center_radius[nearest[0], nearest[1]]
        soma_mask = np.zeros((N, N), dtype=bool)
        rng_soma = np.random.default_rng(SEED + nid*101)
        soma_scale = float(rng_soma.uniform(*CONFIG["soma_scale_range"]))
        soma_angle = float(rng_soma.uniform(0, np.pi))
        sy, sx = ellipse((N//2), (N//2), SOMA_AXES_PX[0]*soma_scale,
                         SOMA_AXES_PX[1]*soma_scale, rotation=soma_angle,
                         shape=(N, N))
        soma_mask[sy, sx] = True
        gt |= soma_mask

        components = ndi.label(gt, structure=np.ones((3,3), dtype=np.uint8))[1]
        if components != 1:
            raise RuntimeError(f"NMO_{nid}: GT contains {components} connected objects.")
        if np.any(gt[0]) or np.any(gt[-1]) or np.any(gt[:,0]) or np.any(gt[:,-1]):
            raise RuntimeError(f"NMO_{nid}: GT touches image boundary.")

        save_tif(OUT / "skeleton" / f"NMO_{nid}_skeleton.tif", skeleton)
        save_tif(OUT / "gt_binary" / f"NMO_{nid}_gt_binary.tif", gt.astype(np.uint8)*255)

        texture = ndi.gaussian_filter(rng_geometry.normal(size=(N,N)), sigma=3)
        texture = texture / max(float(texture.std()), 1e-10)
        shading = ndi.gaussian_filter(rng_geometry.normal(size=(N,N)), sigma=35)
        shading = shading / max(float(shading.std()), 1e-10)
        yy, xx = np.ogrid[:N, :N]
        distance_from_center_um = PX * np.hypot(yy - N//2, xx - N//2)
        proximal_weight = np.exp(-(distance_from_center_um / 10.0)**2)
        pair = []
        for k, (name, p) in enumerate((("clean", CONFIG["clean"]), ("moderate", CONFIG["moderate"]))):
            sigma = p["psf_sigma_px"]
            shape_signal = ndi.gaussian_filter(gt.astype(np.float32), sigma=sigma)
            soma_signal = ndi.gaussian_filter(soma_mask.astype(np.float32), sigma=sigma)
            proximal_signal = ndi.gaussian_filter(
                gt.astype(np.float32) * proximal_weight, sigma=sigma)
            signal = (145 * shape_signal + 90 * soma_signal + 65 * proximal_signal)
            signal *= np.clip(1 + 0.14 * texture, 0.65, 1.35)
            expected = p["background"] + p["bg_variation"]*shading + p["gain"]*signal
            rg = np.random.default_rng(SEED + nid*37 + k)
            noisy = rg.poisson(np.clip(expected, 0, 400)*p["photon_factor"])/p["photon_factor"]
            noisy += rg.normal(0, p["read_noise_sd"], size=(N,N))
            result = np.rint(np.clip(noisy, 0, 255)).astype(np.uint8)
            save_tif(OUT / f"grayscale_{name}" / f"NMO_{nid}_{name}.tif", result)
            pair.append(result)

        qc.append({"neuron_id": nid, "complexity": getattr(row, "complexity_group", "NA"),
                   "source_sha256": sha256(path.read_bytes()).hexdigest(),
                   "width_um": round(float(np.ptp(swc[:,2])),3),
                   "height_um": round(float(np.ptp(swc[:,3])),3),
                   "centerline_margin_px": int(min(rr.min(), cc.min(), N-1-rr.max(), N-1-cc.max())),
                   "width_factor": round(width_factor,4), "soma_scale": round(soma_scale,4),
                   "soma_angle_deg": round(float(np.degrees(soma_angle)),2),
                   "soma_area_gt_um2": round(float(soma_mask.sum())*PX*PX,3),
                   "gt_foreground_pixels": int(gt.sum()),
                   "xy_projected_zero_length_edges": int(sum(np.linalg.norm(xy[i]-xy[j])<1e-8 for j,i in edges)),
                   "gt_connected_components": components})
        preview.append((nid, pair))
        print(f"{len(qc):02d}/30 NMO_{nid}: skeleton + GT binary + Clean + Moderate")

    pd.DataFrame(qc).to_csv(OUT / "qc_30.csv", index=False)
    fig, ax = plt.subplots(6, 10, figsize=(18, 11))
    for idx, (nid, pair) in enumerate(preview):
        for j in (0,1):
            a = ax[idx//5, (idx%5)*2+j]
            a.imshow(pair[j], cmap="gray", vmin=0, vmax=255)
            a.set_title(f"{nid} {'Clean' if j==0 else 'Moderate'}", fontsize=7)
            a.axis("off")
    fig.tight_layout()
    fig.savefig(OUT / "preview_clean_moderate_30.png", dpi=120)
    plt.close(fig)
    print("\nCOMPLETE: 30 skeleton + 30 GT binary + 60 grayscale TIFF")
    print("Folder:", OUT)
    print("Preview:", OUT / "preview_clean_moderate_30.png")
    print("QC:", OUT / "qc_30.csv")
    print("Sholl: measure both GT and Fiji predicted binary in Fiji, identical soma center, 2/2/40 um")
    print("LF morphology v3: same soma ellipse, stronger proximal width and brightness; previous versions preserved")

if __name__ == "__main__":
    main()
