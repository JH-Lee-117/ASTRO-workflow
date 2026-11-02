from pathlib import Path
import argparse as _argparse
_parser = _argparse.ArgumentParser(description="Reproduce the ASTRO SWC preparation step.")
_parser.add_argument('--project-root', type=Path, required=True, help="Analysis project containing data/ and results/.")
_args = _parser.parse_args()
_project = _args.project_root.resolve()
from collections import Counter, defaultdict
from hashlib import sha256
from shutil import copy2
import random
import numpy as np
import pandas as pd

SEED = 20261021
PIXEL_SIZE = 0.151
IMAGE_SIZE = 400
MARGIN_UM = 2.0
HALF_LIMIT = IMAGE_SIZE * PIXEL_SIZE / 2 - MARGIN_UM

SOURCE = _project / "data/swc_all"
DEST = _project / "data/swc_final_low18_medium12"
META = _project / "data/metadata"
MANIFEST = META / "selected_30_low18_medium12.csv"
QC_FILE = META / "all_663_structure_fov_qc.csv"

if MANIFEST.exists():
    existing = pd.read_csv(MANIFEST)
    required = {"neuron_id", "complexity_group", "sha256"}
    if not required.issubset(existing.columns):
        raise RuntimeError("Existing manifest is missing required columns.")
    neuron_ids = pd.to_numeric(existing["neuron_id"], errors="raise")
    if (len(existing) != 30 or not np.isfinite(neuron_ids).all()
            or (neuron_ids <= 0).any()
            or (neuron_ids != np.floor(neuron_ids)).any()
            or neuron_ids.nunique() != 30):
        raise RuntimeError("Existing manifest must contain 30 unique positive integer SWC IDs.")
    groups = existing["complexity_group"]
    if groups.isna().any() or groups.value_counts().to_dict() != {"Low": 18, "Medium": 12}:
        raise RuntimeError("Existing manifest must contain Low=18 and Medium=12.")
    for nid, expected_hash in zip(neuron_ids, existing["sha256"]):
        if (not isinstance(expected_hash, str) or len(expected_hash) != 64
                or any(c not in "0123456789abcdefABCDEF" for c in expected_hash)):
            raise RuntimeError(f"Invalid checksum in manifest for NMO_{int(nid)}.")
        filename = f"NMO_{int(nid)}.swc"
        for folder in (SOURCE, DEST):
            path = folder / filename
            if not path.is_file():
                raise FileNotFoundError(f"Existing selection is incomplete: {path}")
            if sha256(path.read_bytes()).hexdigest() != expected_hash.lower():
                raise RuntimeError(f"Existing selection checksum mismatch: {path}")
    print("Existing final selection verified and preserved:", MANIFEST)
    raise SystemExit(0)

def inspect_swc(path):
    a = np.atleast_2d(np.loadtxt(path, comments="#"))

    if a.shape[1] != 7 or len(a) < 10:
        raise ValueError("Invalid SWC dimensions")
    if not np.isfinite(a).all():
        raise ValueError("Non-finite values")
    if not np.all(a[:, [0, 1, 6]] ==
                  np.round(a[:, [0, 1, 6]])):
        raise ValueError("Invalid node identifiers")

    ids = a[:, 0].astype(int)
    types = a[:, 1].astype(int)
    parents = a[:, 6].astype(int)

    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate node IDs")

    index = {int(n): i for i, n in enumerate(ids)}
    roots = ids[parents == -1]

    if len(roots) != 1:
        raise ValueError("Invalid number of roots")
    if types[index[int(roots[0])]] != 1:
        raise ValueError("Root is not soma")

    children = defaultdict(list)
    for node, parent in zip(ids, parents):
        if parent == -1:
            continue
        if parent not in index:
            raise ValueError("Missing parent")
        children[int(parent)].append(int(node))

    visited = set()
    stack = [int(roots[0])]

    while stack:
        node = stack.pop()
        if node in visited:
            raise ValueError("Cycle detected")
        visited.add(node)
        stack.extend(children[node])

    if len(visited) != len(ids):
        raise ValueError("Disconnected nodes or cycle")

    soma = a[types == 1, 2:4].mean(axis=0)
    xy = a[:, 2:4]
    width, height = np.ptp(xy, axis=0)
    offsets = np.max(np.abs(xy - soma), axis=0)

    bif = sum(
        len(children[int(node)]) >= 2
        for node, kind in zip(ids, types)
        if kind != 1
    )

    return {
        "nodes": len(ids),
        "bifurcations": bif,
        "width_um": float(width),
        "height_um": float(height),
        "z_span_um": float(np.ptp(a[:, 4])),
        "translation_fit": bool(
            width <= 2 * HALF_LIMIT and
            height <= 2 * HALF_LIMIT
        ),
        "soma_center_fit": bool(
            np.all(offsets <= HALF_LIMIT)
        ),
        "sha256": sha256(path.read_bytes()).hexdigest(),
    }

rows = []
files = sorted(SOURCE.glob("NMO_*.swc"))

print("SWC files found:", len(files))

for i, path in enumerate(files, 1):
    nid = int(path.stem.split("_")[1])
    row = {"neuron_id": nid, "status": "PASS"}

    try:
        row.update(inspect_swc(path))
    except Exception as exc:
        row["status"] = "FAIL"
        row["error"] = str(exc)

    rows.append(row)

    if i % 100 == 0:
        print(f"QC checked: {i}/{len(files)}")

qc = pd.DataFrame(rows)
valid = qc[qc.status == "PASS"].copy()
valid = valid.sort_values(
    ["bifurcations", "neuron_id"]
).reset_index(drop=True)

n = len(valid)
valid["complexity_group"] = "High"
valid.loc[:n//3-1, "complexity_group"] = "Low"
valid.loc[n//3:2*n//3-1, "complexity_group"] = "Medium"

qc = qc.merge(
    valid[["neuron_id", "complexity_group"]],
    on="neuron_id",
    how="left"
)
qc.to_csv(QC_FILE, index=False)

print("\nSTRUCTURAL QC")
print("PASS:", len(valid))
print("FAIL:", len(qc) - len(valid))

print("\nFOV QC")
print(
    valid.groupby("complexity_group")[
        ["translation_fit", "soma_center_fit"]
    ].sum().to_string()
)

rng = random.Random(SEED)
selected = []

for group, count in [("Low", 18), ("Medium", 12)]:
    eligible = valid[
        (valid.complexity_group == group) &
        (valid.soma_center_fit)
    ].sort_values("neuron_id")

    if len(eligible) < count:
        raise RuntimeError(
            f"{group}: {len(eligible)} available, "
            f"{count} required"
        )

    indices = rng.sample(list(eligible.index), count)
    selected.append(valid.loc[indices].copy())

final = pd.concat(selected).sort_values(
    "neuron_id"
).reset_index(drop=True)

final["selection_seed"] = SEED
final["pixel_size_um"] = PIXEL_SIZE
final["image_size_px"] = IMAGE_SIZE

DEST.mkdir(parents=True, exist_ok=True)

for nid in final.neuron_id:
    filename = f"NMO_{int(nid)}.swc"
    copy2(SOURCE / filename, DEST / filename)

final.to_csv(MANIFEST, index=False)

print("\nFINAL SELECTION COMPLETE")
print(final.complexity_group.value_counts().to_string())
print("Total:", len(final))
print("FOV:", IMAGE_SIZE * PIXEL_SIZE, "µm square")
print("Manifest:", MANIFEST)
print("SWC folder:", DEST)
