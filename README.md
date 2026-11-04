# ASTRO workflow: GFAP segmentation and Sholl analysis

The **ASTRO workflow** is a semi-automated Fiji/ImageJ procedure for segmenting glial fibrillary acidic protein (GFAP)-positive astrocytic processes and performing soma-centered, image-based Sholl analysis **without manual process tracing**.

This repository provides the **paired-reference validation macro (v6.1.1)** and two Python scripts used to prepare the SWC-derived synthetic GFAP-like images described in the study. These are distinct components:

- **Fiji macro:** segments an 8-bit GFAP grayscale image, runs Sholl analysis, and compares the result with a *supplied* binary reference mask.
- **Python Code1:** quality-checks previously downloaded SWC reconstructions and selects 30 morphologies.
- **Python Code2:** converts the selected reconstructions into 2D skeletons, binary reference masks, and synthetic grayscale images.

**Important:** The current Fiji macro requires a matching binary reference mask for the validation steps. It does not create the reference mask, download data, perform SNT automatic tracing, or run as an unattended batch pipeline. Reference masks are comparison standards, not complete anatomical ground truth.

## Repository files

| File | Purpose |
| --- | --- |
| `ASTRO_workflow_ijm_macro.ijm` | Fiji v6.1.1 paired-reference segmentation and Sholl analysis |
| `Code1_Public_swc_selection_checked.py` | Structural SWC QC, complexity grouping, and reproducible 30-cell selection |
| `Code2_Public_gfap_phantom_renderer_checked.py` | 2D SWC-based binary-reference and grayscale-image generation |
| `rendering_config.json` | Selection seed, spatial calibration, and synthetic-image model settings |
| `LICENSE` | Repository software license |

The Python scripts and configuration operate on a user-supplied project directory; input SWC files and manuscript-level result tables are **not** included in this repository.

## 1. Fiji workflow: paired-reference analysis

### Requirements

Install [Fiji/ImageJ](https://fiji.sc/) with these commands/plugins available:

- **Enhance Local Contrast (CLAHE)**
- **MorphoLibJ** — [Directional Filtering](https://imagej.net/plugins/morpholibj)
- **Auto Threshold** — Triangle and Li thresholding
- **SNT** — [Sholl Analysis (From Image)](https://imagej.net/plugins/snt/)

An inspected Fiji environment used ImageJ **1.54p**, MorphoLibJ **1.6.5**, Auto Threshold **1.18.0**, and SNT **5.0.23**. Compatibility with other versions or operating systems has not been established by this release.

### Prepare one input folder per cell/ROI

```text
my_roi/
├── gfap_image.tif       # Single-channel, single-plane, 8-bit grayscale
└── reference_mask.tif   # Single-plane, 8-bit binary: background=0, foreground=255
```

The reference may instead be a PNG. Both images must have **identical dimensions and spatial alignment**. The macro identifies the reference by its strictly binary (0/255) pixels and the input by its nonbinary grayscale values; it does not identify them from the filenames. Both foreground and background must be present in the binary reference. Keep **exactly two TIFF/PNG image files** in the chosen folder. Input images are not resized or used to generate reference masks.

### Run

1. Open Fiji and save/close any other open image windows.
2. Select **Plugins → Macros → Run...**, then open `ASTRO_workflow_ijm_macro.ijm`.
3. When prompted, select the **folder** containing the grayscale TIFF and matching binary reference.
4. Check the image's spatial calibration. If valid physical calibration is unavailable, supply a verified µm/pixel value or explicitly work in pixel units.
5. Confirm the Sholl settings. The defaults are **start 2, step 2, end 40** in the active units (µm or pixels).
6. Select the **soma center once** on the original grayscale image. The same coordinate is reused for both reference and workflow Sholl analysis.
7. Review the extracted single-cell mask; optionally remove incorrectly connected foreground and complete the interactive QC step.
8. Inspect the saved paired results. If Fiji warns that a Sholl circle extends beyond the image field, treat that circle as **incomplete**, not as a full observation.

### Fiji processing sequence

The grayscale image undergoes background subtraction (50-pixel rolling radius, sliding paraboloid) and CLAHE (block size 127, 256 bins, maximum slope 1.5). Two directional maximum-erosion branches are created (32 orientations):

- **Triangle branch:** 10-pixel directional filter, combined with the CLAHE image at **87.5:12.5** (directional:CLAHE), followed by Triangle thresholding.
- **Li branch:** 5-pixel directional filter, combined at **75:25** (directional:CLAHE), followed by Li thresholding.

The two binary masks are combined using **logical AND**. The macro then retains the largest 8-connected component intersecting a 40-pixel search region around the selected center, followed by visual QC and optional manual removal. Sholl analysis uses the resulting single-cell mask. The **reference mask is not subjected to grayscale preprocessing or segmentation**.

### Fiji output

The macro creates a results directory **inside the selected input folder**:

```text
my_roi/
├── gfap_image.tif
├── reference_mask.tif
└── Results_v6_1_1_BinaryAND/
    ├── Reference/                 # Reference mask and Sholl outputs
    ├── Workflow/                  # Processing intermediates, final mask, Sholl outputs
    ├── Soma_Center.csv
    ├── Soma_Center.roi
    ├── Reference_ShollResults.csv
    ├── Reference_ShollSummary.csv
    ├── Workflow_ShollResults.csv
    ├── Workflow_ShollSummary.csv
    ├── Validation_Comparison.csv
    └── QC_Overlay.tif
```

`Validation_Comparison.csv` includes pixel agreement (Dice, precision, recall), Sholl comparison metrics, calibration, and processing/QC information. In the QC overlay, **white** means overlap, **magenta** workflow-only pixels, **cyan** reference-only pixels, and **black** background.

**No overwrite:** If `Results_v6_1_1_BinaryAND/` already exists, the macro stops rather than replacing results. Rename or archive that directory before another run. The macro generates measurements for **one ROI at a time**, not cohort means or SEM.

## 2. SWC-derived synthetic image preparation (Python)

These scripts reproduce the *selection and image-generation stages* after the source SWCs have been obtained. They **do not** download NeuroMorpho.Org records or reproduce the separate final statistical analyses in the manuscript.

### Python requirements

Use a Python 3 environment with NumPy, pandas, SciPy, scikit-image, tifffile, and Matplotlib. For example, on macOS/Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install numpy pandas scipy scikit-image tifffile matplotlib
```

On Windows, activate the environment with `.venv\Scripts\activate` instead. The `.venv/` directory is local and should not be committed to GitHub.

Package versions are not pinned in this repository; the complete Python pipeline has not been revalidated across arbitrary versions/platforms.

### Obtain the input SWC reconstructions

Source: [NeuroMorpho.Org](https://neuromorpho.org/), **Diniz archive**, GFAP-positive astrocyte reconstructions from mouse hippocampus. The underlying publication is [de Almeida Miranda et al. (2026)](https://doi.org/10.1002/hipo.70085).

For the original study, 666 archive records passed the initial metadata screen. Three pilot reconstructions (`NMO_322485`, `NMO_322659`, `NMO_322774`) were excluded *before* the final structural QC, leaving **663 SWC files**; **660** subsequently passed structural QC. The initial metadata search and downloading were performed separately and are **not implemented by Code1**. To reproduce the reported selection, prepare that same post-pilot set rather than running Code1 on an arbitrary collection of SWCs.

Place downloaded SWC files under `data/swc_all/` and name them `NMO_<neuron_id>.swc` (for example, `NMO_322963.swc`). The scripts expect the original **seven-column SWC data** with physical coordinates in µm. Follow the source archive's data-use and citation requirements.

### Project layout

Run the commands below **from the repository root**, where the Python files and `rendering_config.json` reside:

```text
ASTRO-workflow/
├── README.md
├── LICENSE
├── ASTRO_workflow_ijm_macro.ijm
├── Code1_Public_swc_selection_checked.py
├── Code2_Public_gfap_phantom_renderer_checked.py
├── rendering_config.json
├── data/
│   ├── swc_all/                         # Place the 663 prepared source SWCs here
│   ├── swc_final_low18_medium12/        # Created by Code1
│   └── metadata/                        # Created by Code1
└── results/
    └── phantom_v3_LF/                  # Created by Code2
```

The repository initially contains only the six top-level files; `data/` and `results/` are local working directories created or populated when you run the scripts.

### Step A — SWC structural QC and selection

```bash
python3 Code1_Public_swc_selection_checked.py \
  --project-root . \
  --config rendering_config.json
```

Code1 checks SWC graph structure (node identifiers, soma root, parent links, connectivity, and cycles), ranks structurally valid reconstructions by the number of **non-soma branching nodes**, and divides them into three equal-sized relative-complexity groups. It then selects **18 low- and 12 medium-complexity** reconstructions whose soma-centered projections fit within the image field. The random seed is read from the configuration file (`20261021` in the supplied configuration).

Outputs:

```text
data/metadata/all_663_structure_fov_qc.csv
data/metadata/selected_30_low18_medium12.csv
data/swc_final_low18_medium12/NMO_<ID>.swc
```

The selection manifest includes selected IDs, groups, parameters, and SHA-256 hashes. If the final manifest already exists, Code1 verifies the existing selection and source/copy hashes rather than silently replacing it. **Do not change the candidate SWC collection or settings and assume the same 30 reconstructions will be selected.**

### Step B — Render binary references and GFAP-like grayscale images

```bash
python3 Code2_Public_gfap_phantom_renderer_checked.py \
  --project-root . \
  --config rendering_config.json
```

Code2 verifies the selection manifest and corresponding SWC files. It projects original SWC XY positions onto a **400 × 400-pixel canvas at 0.151 µm/pixel**, preserving physical lengths, and connects the projected nodes using their original SWC parent links. It computes cumulative **3D path distance** along branches and assigns modeled process radii using:

```text
R(d) = f × [1.5 + 9 × exp(−d / 10)] pixels,  f ∈ [0.9, 1.1]
```

It adds a modeled elliptical soma (base semiaxes 24 and 18 pixels; scaling 0.85–1.15), producing an **SWC-derived binary reference**. From each reference, it generates **clean** and **moderate** synthetic GFAP-like grayscale images using spatially varying signals, Gaussian blur, background, and Poisson/Gaussian noise. The moderate images were used for the study's quantitative synthetic-image validation. The process widths, soma dimensions, and signal-generation parameters are modeling choices, **not anatomical measurements of the original cells**.

Outputs:

```text
results/phantom_v3_LF/
├── skeleton/                    # 30 TIFF images
├── gt_binary/                   # 30 binary reference TIFF images
├── grayscale_clean/             # 30 grayscale TIFF images
├── grayscale_moderate/          # 30 grayscale TIFF images
├── qc_30.csv
├── preview_clean_moderate_30.png
└── rendering_config.json        # Configuration snapshot
```

The `phantom_v3_LF` directory name is a historical label for the rendering configuration; **LF is not an experimental condition assigned to the source SWCs**.

The generated `grayscale_moderate/NMO_<ID>_moderate.tif` and corresponding `gt_binary/NMO_<ID>_gt_binary.tif` can be copied **as a pair into a new per-cell folder** for the Fiji macro. Do not put all 60 images in one Fiji input folder.

### Calibration and analysis windows

The provided `rendering_config.json` defines the pixel scale, canvas size, selection seed, and image-model parameters. Some geometry and fluorescence equations are fixed in the Python renderer, so editing a descriptive configuration string **does not change** the implemented equation. Use a separate project/output directory for a changed experimental configuration.

Fiji samples the reference and workflow Sholl profiles at matching radii (typically **2–40 µm in 2-µm steps**). The macro's per-ROI complete-shell MAE/AUC calculations are **not identical** to manuscript-level results computed over fixed common windows:

- **Synthetic Fig. 2:** mean Sholl profiles over **2–28 µm**; presented AUC comparison over **6–28 µm**.
- **Experimental Fig. 3:** Sholl and AUC comparison over **2–16 µm**.

Reproducing group-level manuscript figures therefore requires gathering the per-ROI CSV profiles and applying the stated common-radius windows in a separate analysis step. This repository does not contain that group-summary analysis or the separate **SNT automatic-tracing benchmark**.

## Interpretation and limitations

- The Fiji procedure is **semi-automated**: soma-center selection and mask QC are interactive; a reference mask is required by the supplied *validation* macro.
- The synthetic grayscale images are generated from modeled binary references. This provides a controlled benchmark but does not fully reproduce experimental fluorescence-image complexity.
- The complexity categories are defined **relative to the selected SWC archive pool**, not as universal astrocyte classes. High-complexity reconstructions were excluded from the synthetic benchmark.
- Original SWC spatial geometry is retained, but process thickness and soma dimensions are heuristic modeling parameters.
- Reference masks and ROI-level statistics should not be interpreted as validated whole-cell anatomical boundaries or independent animal-level biological replicates.
- Code1/Code2 reproduce the published **selection and rendering procedures only after the correctly filtered source SWCs have been supplied**. Source acquisition and paper-level statistical aggregation are outside this release.

## Citation and license

If you use the workflow or synthetic preparation scripts, cite the associated ASTRO workflow manuscript when a bibliographic record becomes available, and cite the original SWC data publication and NeuroMorpho.Org as appropriate. The repository software is distributed under the terms in [`LICENSE`](LICENSE); externally sourced images, annotations, and SWC files retain their respective source terms.
