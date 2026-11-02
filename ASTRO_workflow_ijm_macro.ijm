// STAPLE v6.1.1: summary step is corrected only after detailed radius verification; original SNT step is audited.
// ============================================================
// GFAP v6 paired STAPLE reference + DUAL-BRANCH BINARY-AND Sholl analysis
// Fiji / ImageJ Macro Language
// Opens exactly two images from a chosen folder (one gray, one STAPLE binary).
// v6.1 focus/window safeguard 2026-10-23: Triangle[Dir10,CLAHE12.5%] AND Li[Dir5,CLAHE25%]
// Based on user's original v3, keeping grayscale processing and QC intact.
// Accepts any width/height; does not resize input.
// If physical calibration is missing, user must choose pixel units
// or enter a verified micrometers-per-pixel value.
// Retains v5 soma/connected-object/manual-QC/SNT/reference comparison exactly.
// Adds two grayscale-fusion branches, two threshold masks, and binary AND.
// ============================================================


// -------------------- SETTINGS --------------------

// The original 400 x 400 gate has intentionally been removed.
// Images may be 200 x 200, 400 x 400, or other 2D dimensions.
// The source TIFF is never resized by this macro.
// Pixel calibration is determined from TIFF metadata or at runtime.
// The named *_UM parameters below retain original names for compatibility,
// but represent ANALYSIS_UNIT when pixel-only mode is selected.
PIXEL_WIDTH_UM = 1.0;
PIXEL_HEIGHT_UM = 1.0;
ANALYSIS_UNIT = "pixel";
ANALYSIS_UNIT_CSV = "px";

BACKGROUND_RADIUS_PX = 50;
BACKGROUND_SLIDING_PARABOLOID = 1;

CLAHE_BLOCK_SIZE_PX = 127;
CLAHE_HISTOGRAM_BINS = 256;
CLAHE_MAXIMUM_SLOPE = 1.5;
CLAHE_FAST_MODE = 1;

DIRECTIONAL_TYPE = "Max";
DIRECTIONAL_OPERATION = "Erosion";
DIRECTIONAL_DIRECTIONS = 32;
// FIXED v6 research candidate. Change parameters only with a new version.
V6_VERSION = "6.1.1";
V6_FUSION_METHOD = "AND";
DIRECTIONAL_A_LENGTH_PX = 10;
DIRECTIONAL_B_LENGTH_PX = 5;
V6_A_CLAHE_WEIGHT = 0.125;
V6_B_CLAHE_WEIGHT = 0.25;
V6_THRESHOLD_A = "Triangle";
V6_THRESHOLD_B = "Li";
THRESHOLD_IGNORE_BLACK = 1;
THRESHOLD_IGNORE_WHITE = 1;

OBJECT_SEARCH_RADIUS_PX = 40;

SHOLL_START_UM = 2;
SHOLL_STEP_UM  = 2;
SHOLL_END_UM   = 40;
SHOLL_TIMEOUT_MS = 180000;
// Default radii are 2 / 2 / 40 in active analysis units.
// PHANTOM VALIDATION: 400 x 400 at 0.151 um/px; never resize the image.
// They can be changed per run in the setup dialog.
// Keep all compared images on the SAME calibrated radii, center and step.
// For 400px phantoms: radius 40um extends beyond image edges.
// When the retained object has zero pixels on every image edge,
// the extended radius may still be useful for capturing all arbor intersections.
// A warning is shown, and CSV flags partial circles; the user must confirm.
// Fiji/SNT runtime behavior for partial shells still requires empirical QC.

SET_SCALE_SUFFIX = "_setscale.tif";
STEP1_SUFFIX = "_1_BGSubtract_R50px.tif";
STEP2_SUFFIX = "_2_CLAHE_Slope1p5.tif";
STEP3A_SUFFIX = "_3a_Gray_DirectionalMaxErosion10px.tif";
STEP3B_SUFFIX = "_3b_Gray_DirectionalMaxErosion5px.tif";
BLENDA_SUFFIX = "_3c_Fusion_Clahe12p5_Dir10_87p5.tif";
BLENDB_SUFFIX = "_3d_Fusion_Clahe25_Dir5_75.tif";
MASKA_SUFFIX = "_4a_TriangleMask.tif";
MASKB_SUFFIX = "_4b_LiMask.tif";
STEP4_SUFFIX = "_4c_Binary_AND_FinalMask.tif";
SINGLE_SUFFIX = "_5_SingleObject.tif";
ANNOTATED_SUFFIX = "_6_ShollResults.tif";

GEOMETRY_SUFFIX = "_ShollGeometry.csv";
DETAIL_SUFFIX = "_ShollResults.csv";
SUMMARY_SUFFIX = "_ShollSummary.csv";
GRAPH_SUFFIX = "_graph.tif";


// -------------------- PAIRED INPUT / SAFETY --------------------
// Choose a folder holding EXACTLY two images: one non-binary grayscale
// 8-bit TIFF and one 8-bit binary STAPLE TIFF/PNG (0,255 only).
// No input files are rewritten.
// The original v4 closes all image windows on successful completion.
// Safeguard other unrelated Fiji images before opening this pair.
if (nImages > 0)
    exit("Save and close other Fiji image windows before starting this paired macro. It opens the pair itself and closes its images at completion.");
inputDirectory = getDirectory("Select a cell folder (gray TIFF + STAPLE mask)");
inputDirectory = ensureSeparator(inputDirectory);
files = getFileList(inputDirectory);
pairNames = newArray(2);
pairCount = 0;
for (pi=0; pi<files.length; pi++) {
    fn = files[pi];
    low = toLowerCase(fn);
    if (endsWith(low, ".tif") || endsWith(low, ".tiff") ||
        endsWith(low, ".png")) {
        if (pairCount >= 2)
            exit("Folder has more than two image files. Keep exactly one grayscale and one STAPLE mask in this folder.");
        pairNames[pairCount] = fn;
        pairCount++;
    }
}
if (pairCount != 2)
    exit("Expected exactly two image files (one gray TIFF + one binary mask PNG/TIFF). Found: " + pairCount);

grayID = 0;
referenceID = 0;
grayName = "";
referenceName = "";
for (pi=0; pi<2; pi++) {
    open(inputDirectory + pairNames[pi]);
    openedID = getImageID();
    getDimensions(pw, ph, pc, pz, pt);
    if (pc*pz*pt != 1 || bitDepth() != 8)
        exit("Both input images must be single-channel, single-plane, 8-bit. Failed: " + pairNames[pi]);
    getHistogram(hValues, hCounts, 256);
    allPixels = pw*ph;
    isExactBinary = (hCounts[0] + hCounts[255] == allPixels &&
                     hCounts[0] > 0 && hCounts[255] > 0);
    if (isExactBinary) {
        if (referenceID != 0)
            exit("Both images are binary: grayscale TIFF not identified.");
        referenceID = openedID;
        referenceName = pairNames[pi];
        referenceWidth = pw;
        referenceHeight = ph;
    } else {
        if (grayID != 0)
            exit("Neither image is a valid 0/255 binary mask.");
        grayID = openedID;
        grayName = pairNames[pi];
        grayWidth = pw;
        grayHeight = ph;
    }
}
if (grayID == 0 || referenceID == 0)
    exit("Failed to identify the grayscale TIFF and binary reference image.");
if (grayWidth != referenceWidth || grayHeight != referenceHeight)
    exit("Grayscale and STAPLE mask dimensions differ; do not resize here.");
lowGray = toLowerCase(grayName);
if (!(endsWith(lowGray, ".tif") || endsWith(lowGray, ".tiff")))
    exit("The grayscale image must be an 8-bit TIFF.");

resultRoot = inputDirectory + "Results_v6_1_1_BinaryAND" + File.separator;
if (File.exists(resultRoot))
    exit("Results folder already exists: no files were overwritten. Rename or archive it before another run: " + resultRoot);
referenceDirectory = resultRoot + "Reference" + File.separator;
workflowDirectory = resultRoot + "Workflow" + File.separator;

sourceID = grayID;
sourceFileName = grayName;
sourceDirectory = workflowDirectory;
selectImage(sourceID);
getDimensions(imgWidth, imgHeight, imgChannels, imgSlices, imgFrames);
if (imgChannels * imgSlices * imgFrames != 1 || bitDepth() != 8 || is("binary"))
    exit("Grayscale input validation failed.");

if (endsWith(lowGray, ".tiff"))
    BASE = substring(sourceFileName, 0, lengthOf(sourceFileName)-5);
else
    BASE = substring(sourceFileName, 0, lengthOf(sourceFileName)-4);

somaCsvPath = resultRoot + "Soma_Center.csv";
somaRoiPath = resultRoot + "Soma_Center.roi";
referenceResultsPath = resultRoot + "Reference_ShollResults.csv";
referenceSummaryPath = resultRoot + "Reference_ShollSummary.csv";
workflowResultsRootPath = resultRoot + "Workflow_ShollResults.csv";
workflowSummaryRootPath = resultRoot + "Workflow_ShollSummary.csv";
comparisonPath = resultRoot + "Validation_Comparison.csv";
qcOverlayPath = resultRoot + "QC_Overlay.tif";

getPixelSize(
    calUnit,
    currentPixelWidth,
    currentPixelHeight,
    currentVoxelDepth
);

// -------------------- CALIBRATION & SHOLL SETUP --------------------
// Preserve a valid micron calibration instead of force-writing 0.151.
// With unknown calibration, explicitly choose pixel-only or enter um/px.
calibrationNeedsFix = 0;
calibrationSource = "embedded_TIFF";

if (isMicronUnit(calUnit) &&
    currentPixelWidth > 0 &&
    currentPixelHeight > 0) {
    PIXEL_WIDTH_UM = currentPixelWidth;
    PIXEL_HEIGHT_UM = currentPixelHeight;
    ANALYSIS_UNIT = "um";
    ANALYSIS_UNIT_CSV = "um";
} else {
    Dialog.create("Unknown GFAP image calibration");
    Dialog.addMessage(
        "The TIFF does not contain a verified micron scale.\n" +
        "Do not reuse 0.151 um/px for unrelated microscopy data.\n\n" +
        "Image size: " + imgWidth + " x " + imgHeight + " pixels."
    );
    Dialog.addChoice(
        "Analysis scale:",
        newArray("Pixel units (unknown um/px)", "Enter verified um/px"),
        "Pixel units (unknown um/px)"
    );
    Dialog.addNumber("Verified width (um/px), if selected:", 0.151);
    Dialog.addNumber("Verified height (um/px), if selected:", 0.151);
    Dialog.show();
    scaleChoice = Dialog.getChoice();
    scaleWidthInput = Dialog.getNumber();
    scaleHeightInput = Dialog.getNumber();

    calibrationNeedsFix = 1;
    if (scaleChoice == "Enter verified um/px") {
        if (isNaN(scaleWidthInput) || isNaN(scaleHeightInput) ||
            scaleWidthInput <= 0 || scaleHeightInput <= 0)
            exit("Both verified um/px values must be positive.");
        PIXEL_WIDTH_UM = scaleWidthInput;
        PIXEL_HEIGHT_UM = scaleHeightInput;
        ANALYSIS_UNIT = "um";
        ANALYSIS_UNIT_CSV = "um";
        calibrationSource = "user_supplied_um_per_pixel";
    } else {
        PIXEL_WIDTH_UM = 1.0;
        PIXEL_HEIGHT_UM = 1.0;
        ANALYSIS_UNIT = "pixel";
        ANALYSIS_UNIT_CSV = "px";
        calibrationSource = "pixel_only_unverified";
        // Avoid a redundant TIFF duplicate if calibration already is 1 px/px.
        if (toLowerCase(calUnit) == "pixel" &&
            currentPixelWidth == 1 && currentPixelHeight == 1)
            calibrationNeedsFix = 0;
    }
}

Dialog.create("GFAP Sholl settings");
Dialog.addMessage(
    "Image: " + imgWidth + " x " + imgHeight + " px\n" +
    "Scale: " + d2s(PIXEL_WIDTH_UM, 6) + " x " +
    d2s(PIXEL_HEIGHT_UM, 6) + " " + ANALYSIS_UNIT + "/px\n" +
    "Choose identical Sholl radii for every image in a benchmark set."
);
Dialog.addNumber("Starting radius (" + ANALYSIS_UNIT + "):", SHOLL_START_UM);
Dialog.addNumber("Radius step (" + ANALYSIS_UNIT + "):", SHOLL_STEP_UM);
Dialog.addNumber("End radius (" + ANALYSIS_UNIT + "):", SHOLL_END_UM);
Dialog.show();
SHOLL_START_UM = Dialog.getNumber();
SHOLL_STEP_UM = Dialog.getNumber();
SHOLL_END_UM = Dialog.getNumber();

if (isNaN(SHOLL_START_UM) || isNaN(SHOLL_STEP_UM) ||
    isNaN(SHOLL_END_UM) || SHOLL_START_UM < 0 ||
    SHOLL_STEP_UM <= 0 || SHOLL_END_UM <= SHOLL_START_UM)
    exit("Invalid Sholl radius settings.");

intervalCount = (SHOLL_END_UM - SHOLL_START_UM) / SHOLL_STEP_UM;
if (abs(intervalCount - round(intervalCount)) > 0.000001)
    exit("End radius must equal Start + N * Step, for integer N.");

print("=== Input calibration ===");
print("Source: " + sourceFileName);
print("Image size: " + imgWidth + " x " + imgHeight + " px");
print("Calibration source: " + calibrationSource);
print("Pixel scale: " + PIXEL_WIDTH_UM + " x " +
    PIXEL_HEIGHT_UM + " " + ANALYSIS_UNIT + "/pixel");
print("Sholl radii: " + SHOLL_START_UM + " / " +
    SHOLL_STEP_UM + " / " + SHOLL_END_UM +
    " " + ANALYSIS_UNIT);
print("=========================");


// Input folders are left untouched. Require a new Results folder per run.
File.makeDirectory(resultRoot);
File.makeDirectory(referenceDirectory);
File.makeDirectory(workflowDirectory);
if (!File.exists(referenceDirectory) || !File.exists(workflowDirectory))
    exit("Could not create results folders.");

// -------------------- OUTPUT NAMES --------------------

setScaleTitle = BASE + SET_SCALE_SUFFIX;
step1Title = BASE + STEP1_SUFFIX;
step2Title = BASE + STEP2_SUFFIX;
step3ATitle = BASE + STEP3A_SUFFIX;
step3BTitle = BASE + STEP3B_SUFFIX;
blendATitle = BASE + BLENDA_SUFFIX;
blendBTitle = BASE + BLENDB_SUFFIX;
maskATitle = BASE + MASKA_SUFFIX;
maskBTitle = BASE + MASKB_SUFFIX;
step4Title = BASE + STEP4_SUFFIX;
singleTitle = BASE + SINGLE_SUFFIX;
annotatedTitle = BASE + ANNOTATED_SUFFIX;

setScalePath = sourceDirectory + setScaleTitle;
step1Path = sourceDirectory + step1Title;
step2Path = sourceDirectory + step2Title;
step3APath = sourceDirectory + step3ATitle;
step3BPath = sourceDirectory + step3BTitle;
blendAPath = sourceDirectory + blendATitle;
blendBPath = sourceDirectory + blendBTitle;
maskAPath = sourceDirectory + maskATitle;
maskBPath = sourceDirectory + maskBTitle;
step4Path = sourceDirectory + step4Title;
singlePath = sourceDirectory + singleTitle;
annotatedPath = sourceDirectory + annotatedTitle;

geometryPath =
    sourceDirectory + BASE + GEOMETRY_SUFFIX;

resultsPath =
    sourceDirectory + BASE + DETAIL_SUFFIX;

summaryPath =
    sourceDirectory + BASE + SUMMARY_SUFFIX;

graphPath =
    sourceDirectory + BASE + GRAPH_SUFFIX;

scratchTitle =
    BASE + "__component_search";

step3DefaultTitle =
    step2Title + "-directional";

// Temporary filenames generated by Sholl.
tempDetailedTitle =
    singleTitle + "_Sholl-Profiles";

tempDetailedPath =
    sourceDirectory + tempDetailedTitle + ".csv";

tempSummaryPath =
    sourceDirectory + "Sholl_Results.csv";

singleStem =
    substring(
        singleTitle,
        0,
        lengthOf(singleTitle) - 4
    );

tempPlotTitle =
    "Sholl Profile (Linear) for " +
    singleStem + "-1.tif";

tempPlotPath =
    sourceDirectory + tempPlotTitle;


// -------------------- DO NOT OVERWRITE --------------------

if (calibrationNeedsFix)
    abortIfExists(setScalePath);

abortIfExists(step1Path);
abortIfExists(step2Path);
abortIfExists(step3APath);
abortIfExists(step3BPath);
abortIfExists(blendAPath);
abortIfExists(blendBPath);
abortIfExists(maskAPath);
abortIfExists(maskBPath);
abortIfExists(step4Path);
abortIfExists(singlePath);
abortIfExists(geometryPath);
abortIfExists(resultsPath);
abortIfExists(summaryPath);
abortIfExists(graphPath);
abortIfExists(annotatedPath);
abortIfExists(tempDetailedPath);
abortIfExists(tempSummaryPath);
abortIfExists(tempPlotPath);

abortIfOpen(setScaleTitle);
abortIfOpen(step1Title);
abortIfOpen(step2Title);
abortIfOpen(step3ATitle);
abortIfOpen(step3BTitle);
abortIfOpen(blendATitle);
abortIfOpen(blendBTitle);
abortIfOpen(maskATitle);
abortIfOpen(maskBTitle);
abortIfOpen(step4Title);
abortIfOpen(singleTitle);
abortIfOpen(annotatedTitle);
abortIfOpen(scratchTitle);
abortIfOpen(step3DefaultTitle);
abortIfOpen(tempDetailedTitle);
abortIfOpen(tempPlotTitle);


// -------------------- CALIBRATION --------------------

analysisBaseID = sourceID;

if (calibrationNeedsFix) {
    analysisBaseID =
        duplicateFullImage(
            sourceID,
            setScaleTitle
        );

    setVoxelSize(
        PIXEL_WIDTH_UM,
        PIXEL_HEIGHT_UM,
        1.0,
        ANALYSIS_UNIT
    );

    saveTiffOrStop(
        analysisBaseID, setScalePath,
        "calibrated duplicate"
    );
}


// -------------------- SHARED SOMA COORDINATES --------------------
// Point selection is made ONCE, before either Sholl analysis.
selectImage(sourceID);
run("Select None");
setTool("point");
waitForUser("Shared soma center", "Click ONCE at the target soma center in the ORIGINAL GRAYSCALE image.\nThis point will be saved and reused unchanged for the STAPLE direct Sholl and the grayscale workflow.");
selectImage(sourceID);
if (selectionType != 10)
    exit("Exactly one Point ROI is required for the shared soma center.");
getSelectionCoordinates(sharedXs, sharedYs);
if (lengthOf(sharedXs) != 1)
    exit("A single soma point is required.");
somaX = sharedXs[0];
somaY = sharedYs[0];
if (somaX < 0 || somaX >= imgWidth || somaY < 0 || somaY >= imgHeight)
    exit("Shared soma point is outside the image.");
saveAs("Selection", somaRoiPath);
if (!File.exists(somaRoiPath))
    exit("Could not save Soma_Center.roi");
File.saveString("Source_Grayscale,Source_STAPLE,Soma_X_px,Soma_Y_px,Pixel_Width,Pixel_Height,Unit,Start,Step,End\n" +
    grayName + "," + referenceName + "," + d2s(somaX,3) + "," + d2s(somaY,3) + "," +
    d2s(PIXEL_WIDTH_UM,9) + "," + d2s(PIXEL_HEIGHT_UM,9) + "," +
    ANALYSIS_UNIT_CSV + "," + SHOLL_START_UM + "," + SHOLL_STEP_UM + "," + SHOLL_END_UM + "\n", somaCsvPath);
run("Select None");

// -------------------- REFERENCE DIRECT SHOLL --------------------
// Direct reference mode: NO background subtraction, CLAHE, morphology,
// Li threshold, connected-component filtering or manual editing.
refTitle = "STAPLE_Reference.tif";
refPreparedID = duplicateFullImage(referenceID, refTitle);
setVoxelSize(PIXEL_WIDTH_UM, PIXEL_HEIGHT_UM, 1.0, ANALYSIS_UNIT);
refPreparedPath = referenceDirectory + refTitle;
saveTiffOrStop(refPreparedID, refPreparedPath, "unmodified reference binary (calibration metadata only)");
selectImage(refPreparedID);
makePoint(round(somaX), round(somaY));
setSntShollRadiusPreferences(SHOLL_START_UM, SHOLL_STEP_UM, SHOLL_END_UM);
refOptions =
    "datamodechoice=Intersections " +
    "startradius=" + SHOLL_START_UM + " " +
    "stepsize=" + SHOLL_STEP_UM + " " +
    "endradius=" + SHOLL_END_UM + " " +
    "hemishellchoice=[None. Use full shells] " +
    "previewshells=false nspans=1 nspansintchoice=[N/A] " +
    "primarybrancheschoice=[Infer from starting radius] primarybranches=0 " +
    "polynomialchoice=[None. Skip curve fitting] polynomialdegree=0 " +
    "normalizationmethoddescription=[Automatically choose] " +
    "normalizerdescription=Default " +
    "plotoutputdescription=[Linear plot] " +
    "tableoutputdescription=[Detailed & Summary tables] " +
    "annotationsdescription=[ROIs (points and 2D shells)] " +
    "lutchoice=mpl-viridis.lut " +
    "save=true savedir=[" + referenceDirectory + "] " +
    "analysisaction=[Analyze image]";
run("Sholl Analysis (From Image)...", refOptions);
refTempDetailedPath = referenceDirectory + refTitle + "_Sholl-Profiles.csv";
refTempSummaryPath = referenceDirectory + "Sholl_Results.csv";
refTempPlotPath = referenceDirectory + "Sholl Profile (Linear) for STAPLE_Reference-1.tif";
if (waitForThreeFiles(refTempDetailedPath, refTempSummaryPath, refTempPlotPath, SHOLL_TIMEOUT_MS) == 0)
    exit("Reference direct Sholl did not finish or save its CSV/plot files.");
verifyShollRadiiOrStop(refTempDetailedPath, SHOLL_START_UM, SHOLL_STEP_UM, SHOLL_END_UM, 0.000001);
renameFileOrStop(refTempDetailedPath, referenceResultsPath, "Reference detailed Sholl CSV");
compactPairedShollCsv(referenceResultsPath, ANALYSIS_UNIT_CSV);
saveLastCsvRowOrStop(refTempSummaryPath, referenceSummaryPath, "Reference summary Sholl CSV");
renameFileOrStop(refTempPlotPath, referenceDirectory + "Reference_graph.tif", "Reference Sholl graph");

// Re-select original grayscale image for unchanged pre-processing stages.
selectImage(sourceID);
run("Select None");

// -------------------- STEP 1: BACKGROUND SUBTRACTION --------------------

step1ID =
    duplicateFullImage(
        analysisBaseID,
        step1Title
    );

// ROI20 v5 had saved a reference Sholl plot in place of the BG image.
// Select and verify the *created image ID* before AND after the plugin.
selectImage(step1ID);
assertExpectedImageOrStop(step1ID, "before background subtraction");

backgroundOptions =
    "rolling=" + BACKGROUND_RADIUS_PX;

if (BACKGROUND_SLIDING_PARABOLOID)
    backgroundOptions =
        backgroundOptions + " sliding";

// Light background: OFF
// Create background: OFF
// Disable smoothing: OFF
run(
    "Subtract Background...",
    backgroundOptions
);
selectImage(step1ID);
assertExpectedImageOrStop(step1ID, "after background subtraction");

if (bitDepth() != 8)
    exit(
        "Background subtraction failed: " +
        "the output is not 8-bit."
    );

rename(step1Title);

saveTiffOrStop(
    step1ID, step1Path,
    "background-subtracted image"
);


// -------------------- STEP 2: CLAHE --------------------

step2ID =
    duplicateFullImage(
        step1ID,
        step2Title
    );
selectImage(step2ID);
assertExpectedImageOrStop(step2ID, "before CLAHE");

claheOptions =
    "blocksize=" + CLAHE_BLOCK_SIZE_PX +
    " histogram=" + CLAHE_HISTOGRAM_BINS +
    " maximum=" + CLAHE_MAXIMUM_SLOPE +
    " mask=*None*";

if (CLAHE_FAST_MODE)
    claheOptions =
        claheOptions + " fast_(less_accurate)";

run(
    "Enhance Local Contrast (CLAHE)",
    claheOptions
);
selectImage(step2ID);
assertExpectedImageOrStop(step2ID, "after CLAHE");

if (bitDepth() != 8)
    exit(
        "CLAHE failed: " +
        "the output is not 8-bit."
    );

rename(step2Title);

saveTiffOrStop(
    step2ID, step2Path,
    "CLAHE image"
);


// -------------------- STEP 3A/3B: INDEPENDENT DIRECTIONAL GRAYSCALE --------------------
// Both start from the SAME single CLAHE image; no second CLAHE pass.
step3AID = runDirectionalFilter(
    step2ID, step3ATitle, DIRECTIONAL_TYPE, DIRECTIONAL_OPERATION,
    DIRECTIONAL_A_LENGTH_PX, DIRECTIONAL_DIRECTIONS
);
selectImage(step3AID);
if (bitDepth() != 8) exit("Branch A directional grayscale is not 8-bit.");
saveTiffOrStop(step3AID, step3APath, "Branch A directional (10px)");

step3BID = runDirectionalFilter(
    step2ID, step3BTitle, DIRECTIONAL_TYPE, DIRECTIONAL_OPERATION,
    DIRECTIONAL_B_LENGTH_PX, DIRECTIONAL_DIRECTIONS
);
selectImage(step3BID);
if (bitDepth() != 8) exit("Branch B directional grayscale is not 8-bit.");
saveTiffOrStop(step3BID, step3BPath, "Branch B directional (5px)");

// -------------------- STEP 3C/3D: BRANCH-SPECIFIC GRAYSCALE FUSIONS --------------------
// 32-bit weighted addition, display mapped 0..255 before 8-bit conversion.
// Formula A = 0.125 * CLAHE + 0.875 * Directional10.
// Formula B = 0.250 * CLAHE + 0.750 * Directional5.
blendAID = makeWeightedGrayscale(step2ID, step3AID,
    V6_A_CLAHE_WEIGHT, blendATitle);
saveTiffOrStop(blendAID, blendAPath, "Branch A fused grayscale");
blendBID = makeWeightedGrayscale(step2ID, step3BID,
    V6_B_CLAHE_WEIGHT, blendBTitle);
saveTiffOrStop(blendBID, blendBPath, "Branch B fused grayscale");

// -------------------- STEP 4A/4B: SEPARATE AUTO THRESHOLD MASKS --------------------
// Exactly the same Auto Threshold plugin options as v5; different methods.
branchA = makeAutoThresholdMask(blendAID, maskATitle, V6_THRESHOLD_A);
maskAID = branchA[0];
thresholdLowerA = branchA[1];
saveTiffOrStop(maskAID, maskAPath, "Branch A Triangle binary");

branchB = makeAutoThresholdMask(blendBID, maskBTitle, V6_THRESHOLD_B);
maskBID = branchB[0];
thresholdLowerB = branchB[1];
saveTiffOrStop(maskBID, maskBPath, "Branch B Li binary");

// -------------------- STEP 4C: TRUE SPATIAL BINARY AND --------------------
// NOT averaging grayscale cutoffs. NOT averaging the two binary masks.
// Only pixels that are white in BOTH independent masks survive.
imageCalculator("AND create", maskAID, maskBID);
step4ID = getImageID();
if (step4ID == maskAID || step4ID == maskBID)
    exit("Binary AND failed to create a separate output window.");
selectImage(step4ID);
rename(step4Title);
setVoxelSize(PIXEL_WIDTH_UM, PIXEL_HEIGHT_UM, 1.0, ANALYSIS_UNIT);
binaryQcOrStop("Final binary AND mask");
saveTiffOrStop(step4ID, step4Path, "Final v6 binary AND mask");

print("v6 branch A: Dir" + DIRECTIONAL_A_LENGTH_PX +
    " CLAHE=" + V6_A_CLAHE_WEIGHT + " " + V6_THRESHOLD_A +
    " cutoff=" + thresholdLowerA);
print("v6 branch B: Dir" + DIRECTIONAL_B_LENGTH_PX +
    " CLAHE=" + V6_B_CLAHE_WEIGHT + " " + V6_THRESHOLD_B +
    " cutoff=" + thresholdLowerB);
print("v6 binary combination: " + V6_FUSION_METHOD);


// -------------------- STEP 5: SAME SAVED SOMA CENTER --------------------
// somaX/somaY were selected once and saved BEFORE reference direct Sholl.
// Do not solicit a second point; keep both analyses in identical coordinates.
selectImage(step4ID);
run("Select None");

// -------------------- STEP 6: SINGLE OBJECT --------------------

largest = findLargestComponentNear(
    step4ID,
    somaX,
    somaY,
    imgWidth,
    imgHeight,
    OBJECT_SEARCH_RADIUS_PX,
    scratchTitle
);

if (largest[3] == 0)
    exit(
        "No white component intersects " +
        "the " + OBJECT_SEARCH_RADIUS_PX + " px search circle."
    );

seedX = largest[0];
seedY = largest[1];

singleID =
    duplicateFullImage(
        step4ID,
        singleTitle
    );

objectPixels =
    retainWhiteComponent(
        seedX,
        seedY
    );

rename(singleTitle);

if (objectPixels <= 0)
    exit(
        "The retained object contains no white pixels."
    );

initialSingleObjectPixels =
    objectPixels;


// -------------------- STEP 7: QC, EDIT & RETRY --------------------

qcAccepted = 0;
manualEditRounds = 0;
manualPixelsRemoved = 0;
disconnectedPixelsRemoved = 0;

while (qcAccepted == 0) {
    qcDecision =
        askSingleObjectQc();

    if (startsWith(qcDecision, "PASS")) {
        qcAccepted = 1;
    } else {
        removedThisRound =
            manuallyEditSingleObject(
                singleID,
                singleTitle
            );

        manualPixelsRemoved =
            manualPixelsRemoved +
            removedThisRound;

        manualEditRounds++;

        selectImage(singleID);
        binaryQcOrStop(singleTitle);

        beforeDebrisCleanup =
            countWhitePixels();

        cleaned =
            retainSomaConnectedComponent(
                singleID,
                somaX,
                somaY,
                seedX,
                seedY,
                imgWidth,
                imgHeight,
                OBJECT_SEARCH_RADIUS_PX
            );

        seedX = cleaned[0];
        seedY = cleaned[1];
        objectPixels = cleaned[2];

        removedDisconnectedThisRound =
            beforeDebrisCleanup -
            objectPixels;

        disconnectedPixelsRemoved =
            disconnectedPixelsRemoved +
            removedDisconnectedThisRound;

        binaryQcOrStop(singleTitle);
    }
}


// -------------------- SAVE SINGLE OBJECT --------------------

if (manualEditRounds > 0)
    manualEditPerformed = 1;
else
    manualEditPerformed = 0;

totalPixelsRemoved =
    initialSingleObjectPixels -
    objectPixels;

if (totalPixelsRemoved < 0)
    exit(
        "The edited SingleObject contains " +
        "additional white pixels."
    );

if (
    manualPixelsRemoved +
    disconnectedPixelsRemoved !=
    totalPixelsRemoved
)
    exit(
        "Manual/debris pixel accounting did not match."
    );

selectImage(singleID);
run("Select None");

saveTiffOrStop(
    singleID, singlePath,
    "single-object binary mask"
);


// -------------------- GEOMETRY --------------------

imageGeometry = getImageGeometry(
    somaX,
    somaY,
    imgWidth,
    imgHeight,
    PIXEL_WIDTH_UM,
    PIXEL_HEIGHT_UM
);

largestCompleteShell =
    imageGeometry[0];

farthestImageCorner =
    imageGeometry[1];

objectMaxRadius = getObjectMaxRadius(
    somaX,
    somaY,
    imgWidth,
    imgHeight,
    PIXEL_WIDTH_UM,
    PIXEL_HEIGHT_UM
);

objectTouchesEdge =
    whiteTouchesImageEdge(
        imgWidth,
        imgHeight
    );

// Do not silently treat incomplete concentric shells as full Sholl data.
shollCompleteShell = 1;
if (SHOLL_END_UM > largestCompleteShell + 0.000001) {
    shollCompleteShell = 0;
    Dialog.create("WARNING: incomplete Sholl circles");
    Dialog.addMessage(
        "Requested end radius: " + SHOLL_END_UM + " " + ANALYSIS_UNIT + "\n" +
        "Largest COMPLETE shell: " + d2s(largestCompleteShell, 3) +
        " " + ANALYSIS_UNIT + "\n\n" +
        "Circles beyond this limit extend outside the image.\n" +
        "The full circle is not observed: use the CSV partial-shell flag.\n" +
        "Confirm that the target object does NOT touch the image edge."
    );
    Dialog.addChoice(
        "Action:",
        newArray("STOP and choose a smaller radius", "PROCEED exploratory (flag CSV)"),
        "PROCEED exploratory (flag CSV)"
    );
    Dialog.show();
    shellChoice = Dialog.getChoice();
    if (shellChoice == "STOP and choose a smaller radius") {
        // These intermediate outputs were created by this run only;
        // every path was checked for preexistence before the first write.
        if (File.exists(step1Path)) File.delete(step1Path);
        if (File.exists(step2Path)) File.delete(step2Path);
        if (File.exists(step3APath)) File.delete(step3APath);
        if (File.exists(step3BPath)) File.delete(step3BPath);
        if (File.exists(blendAPath)) File.delete(blendAPath);
        if (File.exists(blendBPath)) File.delete(blendBPath);
        if (File.exists(maskAPath)) File.delete(maskAPath);
        if (File.exists(maskBPath)) File.delete(maskBPath);
        if (File.exists(step4Path)) File.delete(step4Path);
        if (File.exists(singlePath)) File.delete(singlePath);
        if (calibrationNeedsFix && File.exists(setScalePath))
            File.delete(setScalePath);
        exit("Stopped: end radius exceeds largest complete Sholl shell.\n" +
             "The intermediate TIFF outputs were removed; source unchanged.");
    }
}

print(
    "=== Sholl geometry: " +
    sourceFileName + " ==="
);

print(
    "Largest complete shell: " +
    d2s(largestCompleteShell, 3) +
    " " + ANALYSIS_UNIT
);

print(
    "Farthest image corner: " +
    d2s(farthestImageCorner, 3) +
    " " + ANALYSIS_UNIT
);

print(
    "Retained-object max radius: " +
    d2s(objectMaxRadius, 3) +
    " " + ANALYSIS_UNIT
);

print(
    "Object touches edge: " +
    objectTouchesEdge
);

print(
    "Manual edit rounds: " +
    manualEditRounds
);

print(
    "Manual pixels removed: " +
    manualPixelsRemoved
);

print(
    "Disconnected pixels removed: " +
    disconnectedPixelsRemoved
);

print(
    "Total pixels removed: " +
    totalPixelsRemoved
);

print("========================================");


// -------------------- GEOMETRY CSV --------------------

geometryHeader =
    "Source_File," +
    "Soma_X_px," +
    "Soma_Y_px," +
    "Largest_Complete_Shell_" + ANALYSIS_UNIT_CSV + "," +
    "Farthest_Image_Corner_" + ANALYSIS_UNIT_CSV + "," +
    "Object_Max_Radius_" + ANALYSIS_UNIT_CSV + "," +
    "Object_Touches_Edge," +
    "Image_Width_px,Image_Height_px," +
    "Pixel_Width_in_Analysis_Units,Pixel_Height_in_Analysis_Units," +
    "Analysis_Unit,Calibration_Source,Sholl_End_Within_Complete_Shell," +
    "Retained_Object_Pixels," +
    "CLAHE_Block_Size_px," +
    "CLAHE_Histogram_Bins," +
    "CLAHE_Maximum_Slope," +
    "v6_Version,v6_Binary_Fusion," +
    "BranchA_Directional_px,BranchA_CLAHE_Weight,BranchA_Directional_Weight," +
    "BranchA_Threshold_Method,BranchA_Lower_Threshold," +
    "BranchB_Directional_px,BranchB_CLAHE_Weight,BranchB_Directional_Weight," +
    "BranchB_Threshold_Method,BranchB_Lower_Threshold," +
    "Background_Radius_px," +
    "Threshold_Method," +
    "Ignore_Black," +
    "Ignore_White," +
    "Manual_Edit_Performed," +
    "Manual_Edit_Rounds," +
    "Initial_SingleObject_Pixels," +
    "Manual_Pixels_Removed," +
    "Disconnected_Pixels_Removed," +
    "Total_Pixels_Removed," +
    "Sholl_Start_" + ANALYSIS_UNIT_CSV + "," +
    "Sholl_Step_" + ANALYSIS_UNIT_CSV + "," +
    "Sholl_End_" + ANALYSIS_UNIT_CSV + "\n";

quotedSourceFileName =
    replace(
        sourceFileName,
        "\"",
        "\"\""
    );

quotedSourceFileName =
    "\"" +
    quotedSourceFileName +
    "\"";

geometryRow =
    quotedSourceFileName + "," +
    d2s(somaX, 3) + "," +
    d2s(somaY, 3) + "," +
    d2s(largestCompleteShell, 6) + "," +
    d2s(farthestImageCorner, 6) + "," +
    d2s(objectMaxRadius, 6) + "," +
    objectTouchesEdge + "," +
    imgWidth + "," + imgHeight + "," +
    d2s(PIXEL_WIDTH_UM, 9) + "," +
    d2s(PIXEL_HEIGHT_UM, 9) + "," +
    ANALYSIS_UNIT + "," + calibrationSource + "," +
    shollCompleteShell + "," +
    objectPixels + "," +
    CLAHE_BLOCK_SIZE_PX + "," +
    CLAHE_HISTOGRAM_BINS + "," +
    CLAHE_MAXIMUM_SLOPE + "," +
    V6_VERSION + "," + V6_FUSION_METHOD + "," +
    DIRECTIONAL_A_LENGTH_PX + "," + V6_A_CLAHE_WEIGHT + "," +
    (1.0-V6_A_CLAHE_WEIGHT) + "," + V6_THRESHOLD_A + "," + thresholdLowerA + "," +
    DIRECTIONAL_B_LENGTH_PX + "," + V6_B_CLAHE_WEIGHT + "," +
    (1.0-V6_B_CLAHE_WEIGHT) + "," + V6_THRESHOLD_B + "," + thresholdLowerB + "," +
    BACKGROUND_RADIUS_PX + "," +
    "Triangle_AND_Li" + "," +
    THRESHOLD_IGNORE_BLACK + "," +
    THRESHOLD_IGNORE_WHITE + "," +
    manualEditPerformed + "," +
    manualEditRounds + "," +
    initialSingleObjectPixels + "," +
    manualPixelsRemoved + "," +
    disconnectedPixelsRemoved + "," +
    totalPixelsRemoved + "," +
    SHOLL_START_UM + "," +
    SHOLL_STEP_UM + "," +
    SHOLL_END_UM + "\n";

File.saveString(
    geometryHeader +
    geometryRow,
    geometryPath
);

if (File.exists(geometryPath) == 0)
    exit(
        "Saving failed for the geometry CSV."
    );


// -------------------- STEP 8: SHOLL --------------------

selectImage(singleID);

makePoint(
    round(somaX),
    round(somaY)
);

setSntShollRadiusPreferences(
    SHOLL_START_UM,
    SHOLL_STEP_UM,
    SHOLL_END_UM
);

shollOptions =
    "datamodechoice=Intersections " +
    "startradius=" + SHOLL_START_UM + " " +
    "stepsize=" + SHOLL_STEP_UM + " " +
    "endradius=" + SHOLL_END_UM + " " +
    "hemishellchoice=[None. Use full shells] " +
    "previewshells=false " +
    "nspans=1 " +
    "nspansintchoice=[N/A] " +
    "primarybrancheschoice=[Infer from starting radius] " +
    "primarybranches=0 " +
    "polynomialchoice=[None. Skip curve fitting] " +
    "polynomialdegree=0 " +
    "normalizationmethoddescription=[Automatically choose] " +
    "normalizerdescription=Default " +
    "plotoutputdescription=[Linear plot] " +
    "tableoutputdescription=[Detailed & Summary tables] " +
    "annotationsdescription=[ROIs (points and 2D shells)] " +
    "lutchoice=mpl-viridis.lut " +
    "save=true " +
    "savedir=[" + sourceDirectory + "] " +
    "analysisaction=[Analyze image]";

run(
    "Sholl Analysis (From Image)...",
    shollOptions
);

if (
    waitForThreeFiles(
        tempDetailedPath,
        tempSummaryPath,
        tempPlotPath,
        SHOLL_TIMEOUT_MS
    ) == 0
)
    exit(
        "Sholl analysis or saving did not finish " +
        "within the timeout."
    );

verifyShollRadiiOrStop(
    tempDetailedPath,
    SHOLL_START_UM,
    SHOLL_STEP_UM,
    SHOLL_END_UM,
    0.000001
);

renameFileOrStop(
    tempDetailedPath,
    resultsPath,
    "Detailed Sholl CSV"
);

// -------------------- COMPACT SHOLL CSV --------------------

// Read the saved CSV itself. No Results-window title is assumed.
compactText = replace(
    File.openAsString(resultsPath),
    "\r",
    ""
);

compactLines = split(compactText, "\n");
compactRadii = newArray(compactLines.length);
compactIntersections = newArray(compactLines.length);

compactRowCount = 0;
shollAUC = 0;
shollMaxIntersections = -1;
shollEnclosingRadius = 0;

for (
    compactLine = 1;
    compactLine < compactLines.length;
    compactLine++
) {
    if (lengthOf(compactLines[compactLine]) > 0) {
        compactColumns = split(
            compactLines[compactLine],
            ","
        );

        // Sholl CSV columns: row label, Radius, Inters., ...
        if (compactColumns.length >= 3) {
            compactRadius = parseFloat(compactColumns[1]);
            compactIntersection = parseFloat(compactColumns[2]);

            compactRadii[compactRowCount] = compactRadius;
            compactIntersections[compactRowCount] =
                compactIntersection;

            if (compactIntersection > shollMaxIntersections)
                shollMaxIntersections = compactIntersection;

            if (compactIntersection > 0)
                shollEnclosingRadius = compactRadius;

            if (compactRowCount > 0) {
                compactRadiusInterval =
                    compactRadius -
                    compactRadii[compactRowCount - 1];

                shollAUC =
                    shollAUC +
                    (
                        compactIntersections[compactRowCount - 1] +
                        compactIntersection
                    ) / 2 * compactRadiusInterval;
            }

            compactRowCount++;
        }
    }
}

if (compactRowCount < 2)
    exit(
        "The saved Sholl CSV does not contain enough rows."
    );

compactOutput =
    "Radius_" + ANALYSIS_UNIT_CSV +
    ",Intersections,Summary_Metric,Summary_Value\n";

for (
    compactRow = 0;
    compactRow < compactRowCount;
    compactRow++
) {
    compactMetric = "";
    compactValue = "";

    if (compactRow == 0) {
        compactMetric = "Sholl_AUC";
        compactValue = d2s(shollAUC, 3);
    }

    if (compactRow == 1) {
        compactMetric = "Max_Intersections";
        compactValue = d2s(shollMaxIntersections, 0);
    }

    if (compactRow == 2) {
        compactMetric = "Enclosing_Radius_" + ANALYSIS_UNIT_CSV;
        compactValue = d2s(shollEnclosingRadius, 3);
    }

    compactOutput =
        compactOutput +
        d2s(compactRadii[compactRow], 3) + "," +
        d2s(compactIntersections[compactRow], 0) + "," +
        compactMetric + "," +
        compactValue + "\n";
}

File.saveString(compactOutput, resultsPath);

compactCheck = File.openAsString(resultsPath);

if (indexOf(compactCheck, "Radius_" + ANALYSIS_UNIT_CSV +
    ",Intersections") != 0)
    exit(
        "Saving the compact Sholl CSV failed."
    );

saveLastCsvRowOrStop(
    tempSummaryPath,
    summaryPath,
    "Summary Sholl CSV"
);

renameFileOrStop(
    tempPlotPath,
    graphPath,
    "Sholl graph TIFF"
);


// -------------------- STEP 9: ANNOTATED IMAGE --------------------

selectImage(singleID);

if (Overlay.size() == 0)
    exit(
        "No Sholl annotation overlay was found."
    );

run("Select None");

imageCountBeforeFlatten =
    nImages;

run("Flatten");

if (nImages != imageCountBeforeFlatten + 1)
    exit(
        "Could not create the flattened " +
        "annotation image."
    );

annotatedID = getImageID();
rename(annotatedTitle);

saveTiffOrStop(
    annotatedID, annotatedPath,
    "flattened Sholl annotation image"
);


// -------------------- PAIRED VALIDATION --------------------
// Preserve the existing full workflow outputs and also provide a flat
// pair of result tables and a colored per-pixel agreement QC image.
if (File.exists(workflowResultsRootPath) || File.exists(workflowSummaryRootPath))
    exit("Paired workflow result already exists: overwrite forbidden.");
File.copy(resultsPath, workflowResultsRootPath);
File.copy(summaryPath, workflowSummaryRootPath);
if (!File.exists(workflowResultsRootPath) || !File.exists(workflowSummaryRootPath))
    exit("Copying paired workflow Sholl CSVs failed.");

// Read both masks (0/255) into arrays WITHOUT modifying either image.
refBits = newArray(imgWidth*imgHeight);
workBits = newArray(imgWidth*imgHeight);
selectImage(refPreparedID);
for (qy=0; qy<imgHeight; qy++) {
    for (qx=0; qx<imgWidth; qx++)
        if (getPixel(qx,qy) == 255) refBits[qy*imgWidth+qx] = 1;
        else refBits[qy*imgWidth+qx] = 0;
}
selectImage(singleID);
for (qy=0; qy<imgHeight; qy++) {
    for (qx=0; qx<imgWidth; qx++)
        if (getPixel(qx,qy) == 255) workBits[qy*imgWidth+qx] = 1;
        else workBits[qy*imgWidth+qx] = 0;
}

// TIFF color key: black TN, white TP, MAGENTA FP, CYAN FN.
newImage("QC_Overlay", "RGB black", imgWidth, imgHeight, 1);
qcID = getImageID();
TP=0; FP=0; FN=0; TN=0;
for (qy=0; qy<imgHeight; qy++) {
    for (qx=0; qx<imgWidth; qx++) {
        qidx = qy*imgWidth+qx;
        rbit=refBits[qidx];
        wbit=workBits[qidx];
        if (rbit==1 && wbit==1) {TP++; setPixel(qx,qy,16777215);}
        else if (rbit==0 && wbit==1) {FP++; setPixel(qx,qy,16711935);}
        else if (rbit==1 && wbit==0) {FN++; setPixel(qx,qy,65535);}
        else {TN++;}
    }
}
saveTiffOrStop(qcID, qcOverlayPath, "RGB QC_Overlay");

precision=0; recall=0; dice=0;
if (TP+FP>0) precision=TP/(TP+FP);
if (TP+FN>0) recall=TP/(TP+FN);
if (2*TP+FP+FN>0) dice=2*TP/(2*TP+FP+FN);

// Compare identical radii, using only COMPLETE shells for primary errors.
refText=replace(File.openAsString(referenceResultsPath),"\r","");
workText=replace(File.openAsString(workflowResultsRootPath),"\r","");
refLines=split(refText,"\n");
workLines=split(workText,"\n");
if (refLines.length!=workLines.length)
    exit("Paired Sholl CSV row counts do not match.");
validN=0; absSum=0; sqSum=0; validRefAUC=0; validWorkAUC=0;
prevR=0; prevRef=0; prevWork=0;
validMaxRadius = largestCompleteShell;
for (mi=1;mi<refLines.length;mi++) {
    if (lengthOf(refLines[mi])>0) {
        if (lengthOf(workLines[mi])==0) exit("Paired Sholl CSV is missing a row.");
        rc=split(refLines[mi],","); wc=split(workLines[mi],",");
        if (rc.length<2 || wc.length<2) exit("Malformed paired Sholl CSV.");
        rr=parseFloat(rc[0]); wr=parseFloat(wc[0]);
        rv=parseFloat(rc[1]); wv=parseFloat(wc[1]);
        if (abs(rr-wr)>0.000001) exit("Paired Sholl radii do not match.");
        if (rr<=validMaxRadius+0.000001) {
            validN++;
            absSum+=abs(rv-wv);
            sqSum+=(rv-wv)*(rv-wv);
            if (validN>1) {
                validRefAUC+=(prevRef+rv)/2*(rr-prevR);
                validWorkAUC+=(prevWork+wv)/2*(rr-prevR);
            }
            prevR=rr;prevRef=rv;prevWork=wv;
        }
    }
}
mae=-1; rmse=-1; aucErrorPct=-1;
if(validN>0){mae=absSum/validN;rmse=sqrt(sqSum/validN);}
if(validN>1 && validRefAUC>0)
    aucErrorPct=100*abs(validWorkAUC-validRefAUC)/validRefAUC;

cmpHeader="Gray_File,Reference_File,Soma_X_px,Soma_Y_px,Unit,Start,Step,End,Largest_Complete_Shell,Complete_Radii_N,Sholl_MAE_Complete,Sholl_RMSE_Complete,Reference_AUC_Complete,Workflow_AUC_Complete,AUC_Absolute_Error_Percent,TP,FP,FN,TN,Dice,Precision,Recall,Manual_Edit_Rounds,Manual_Pixels_Removed,Object_Touches_Edge,Color_Key,v6_Version,v6_Mask_Fusion,BranchA_Threshold,BranchB_Threshold,DirA_px,DirB_px,CLAHE_A_Weight,CLAHE_B_Weight,Threshold_A_Lower,Threshold_B_Lower\n";
cmpRow=grayName+","+referenceName+","+d2s(somaX,3)+","+d2s(somaY,3)+","+ANALYSIS_UNIT_CSV+","+
    SHOLL_START_UM+","+SHOLL_STEP_UM+","+SHOLL_END_UM+","+d2s(validMaxRadius,6)+","+validN+","+
    d2s(mae,6)+","+d2s(rmse,6)+","+d2s(validRefAUC,6)+","+d2s(validWorkAUC,6)+","+
    d2s(aucErrorPct,6)+","+TP+","+FP+","+FN+","+TN+","+d2s(dice,6)+","+d2s(precision,6)+","+
    d2s(recall,6)+","+manualEditRounds+","+manualPixelsRemoved+","+objectTouchesEdge+",TP_white_FP_magenta_FN_cyan_TN_black,"+
    V6_VERSION+","+V6_FUSION_METHOD+","+V6_THRESHOLD_A+","+V6_THRESHOLD_B+","+
    DIRECTIONAL_A_LENGTH_PX+","+DIRECTIONAL_B_LENGTH_PX+","+
    V6_A_CLAHE_WEIGHT+","+V6_B_CLAHE_WEIGHT+","+
    thresholdLowerA+","+thresholdLowerB+"\n";
File.saveString(cmpHeader+cmpRow,comparisonPath);
if (!File.exists(comparisonPath)) exit("Could not save Validation_Comparison.csv");

// -------------------- COMPLETE --------------------

createdFiles =
    step1Title + "\n" +
    step2Title + "\n" +
    step3ATitle + "\n" +
    step3BTitle + "\n" +
    blendATitle + "\n" +
    blendBTitle + "\n" +
    maskATitle + "\n" +
    maskBTitle + "\n" +
    step4Title + "\n" +
    singleTitle + "\n" +
    BASE + GEOMETRY_SUFFIX + "\n" +
    BASE + DETAIL_SUFFIX + "\n" +
    BASE + SUMMARY_SUFFIX + "\n" +
    BASE + GRAPH_SUFFIX + "\n" +
    annotatedTitle;

if (calibrationNeedsFix)
    createdFiles =
        setScaleTitle + "\n" +
        createdFiles;

showStatus(
    "Single-cell Sholl analysis completed: QC PASS"
);

showMessage(
    "Single-cell Sholl analysis complete",
    "Original file: " +
    sourceFileName + "\n" +
    "Soma center: " +
    d2s(somaX, 2) + ", " +
    d2s(somaY, 2) + " px\n" +
    "Selected seed: " +
    seedX + ", " +
    seedY + " px\n" +
    "Retained object: " +
    objectPixels + " pixels\n" +
    "Manual edit rounds: " +
    manualEditRounds + "\n" +
    "Manual pixels removed: " +
    manualPixelsRemoved + "\n" +
    "Disconnected pixels removed: " +
    disconnectedPixelsRemoved + "\n" +
    "Object max radius: " +
    d2s(objectMaxRadius, 2) + " " + ANALYSIS_UNIT + "\n" +
    "Largest complete shell: " +
    d2s(largestCompleteShell, 2) + " " + ANALYSIS_UNIT + "\n" +
    "Sholl AUC: " +
    d2s(shollAUC, 3) + " intersections*" + ANALYSIS_UNIT + "\n" +
    "Maximum intersections: " +
    d2s(shollMaxIntersections, 0) + "\n" +
    "Enclosing radius: " +
    d2s(shollEnclosingRadius, 2) + " " + ANALYSIS_UNIT + "\n" +
    "Sholl full-shell QC: " + shollCompleteShell + "\n" +
    "Sholl radii: " +
    SHOLL_START_UM + " / " +
    SHOLL_STEP_UM + " / " +
    SHOLL_END_UM + " " + ANALYSIS_UNIT + "\n" +
    "Output folder: " +
    sourceDirectory + "\n\n" +
    "Created files:\n" +
    createdFiles + "\n\n" +
    "QC status: PASS"
);


// -------------------- CLOSE WINDOWS --------------------

// All required outputs have been saved successfully.
run("Close All");

closeTitles = getList("window.titles");

for (
    closeIndex = closeTitles.length - 1;
    closeIndex >= 0;
    closeIndex--
) {
    closeTitle = closeTitles[closeIndex];
    closeLower = toLowerCase(closeTitle);

    closeThisWindow = 0;

    if (closeLower == "results")
        closeThisWindow = 1;

    if (closeLower == "log")
        closeThisWindow = 1;

    if (indexOf(closeLower, "sholl") >= 0)
        closeThisWindow = 1;

    if (indexOf(closeLower, "roi manager") >= 0)
        closeThisWindow = 1;

    if (indexOf(closeLower, ".ijm") >= 0)
        closeThisWindow = 0;

    if (closeThisWindow) {
        selectWindow(closeTitle);
        run("Close");
    }
}


// -------------------- PHANTOM VALIDATION: NO SOMA AREA --------------------
// Original v4 Soma-area opening/measurement intentionally removed.
// This avoids the post-Sholl binary-QC failure when opening radius 12 px
// eliminates small phantom somata. Soma-clicking, selection, Sholl,
// geometry CSV and image/csv exports remain unchanged.

// ============================================================
// HELPER FUNCTIONS
// ============================================================

// Compact reference using the SAME metric definitions as the v4 workflow.
function compactPairedShollCsv(csvPath, unitCsv) {
    raw = replace(File.openAsString(csvPath), "\r", "");
    lines = split(raw, "\n");
    radii=newArray(lines.length);
    intersections=newArray(lines.length);
    n=0; auc=0; maxN=-1; enclosing=0;
    for (i=1; i<lines.length; i++) {
        if (lengthOf(lines[i])>0) {
            columns=split(lines[i],",");
            if (columns.length>=3) {
                radius=parseFloat(columns[1]);
                hits=parseFloat(columns[2]);
                radii[n]=radius; intersections[n]=hits;
                if (hits>maxN) maxN=hits;
                if (hits>0) enclosing=radius;
                if (n>0) auc+=(intersections[n-1]+hits)/2*(radius-radii[n-1]);
                n++;
            }
        }
    }
    if (n<2) exit("Reference Sholl returned too few radius rows.");
    out="Radius_"+unitCsv+",Intersections,Summary_Metric,Summary_Value\n";
    for(i=0;i<n;i++) {
        metric=""; val="";
        if(i==0){metric="Sholl_AUC";val=d2s(auc,3);}
        if(i==1){metric="Max_Intersections";val=d2s(maxN,0);}
        if(i==2){metric="Enclosing_Radius_"+unitCsv;val=d2s(enclosing,3);}
        out+=d2s(radii[i],3)+","+d2s(intersections[i],0)+","+metric+","+val+"\n";
    }
    File.saveString(out,csvPath);
}

function ensureSeparator(path) {
    if (endsWith(path, File.separator) == 0)
        path = path + File.separator;

    return path;
}

function isMicronUnit(unit) {
    unit = toLowerCase(unit);

    if (unit == "µm") return 1;
    if (unit == "μm") return 1;
    if (unit == "um") return 1;
    if (unit == "micron") return 1;
    if (unit == "microns") return 1;
    if (unit == "micrometer") return 1;
    if (unit == "micrometers") return 1;

    return 0;
}

function abortIfExists(path) {
    if (File.exists(path))
        exit(
            "Do not overwrite: output already exists.\n\n" +
            path
        );
}

function abortIfOpen(title) {
    if (isOpen(title))
        exit(
            "A required output window is already open.\n\n" +
            title
        );
}

function duplicateFullImage(
    sourceImageID,
    newTitle
) {
    selectImage(sourceImageID);

    hadSelection =
        selectionType != -1;

    if (hadSelection)
        Roi.copy();

    getDimensions(w, h, c, z, t);
    if (w != imgWidth || h != imgHeight || c*z*t != 1)
        exit("Duplicate source geometry mismatch for " + newTitle +
             ". Actual " + w + "x" + h + " C/Z/T=" + c + "/" + z + "/" + t +
             "; expected " + imgWidth + "x" + imgHeight + " single plane.");

    makeRectangle(
        0,
        0,
        w,
        h
    );

    run(
        "Duplicate...",
        "title=[" + newTitle + "]"
    );

    // SNT may bring a graph window to the front. Do NOT infer the new ID
    // from whichever window happens to be active immediately after Duplicate.
    if (!isOpen(newTitle))
        exit("Duplicate window not found: " + newTitle);
    selectWindow(newTitle);
    duplicateID = getImageID();
    if (duplicateID == sourceImageID)
        exit("Duplicate reused the source image ID: " + newTitle);
    assertExpectedImageOrStop(duplicateID, "newly duplicated " + newTitle);

    selectImage(sourceImageID);

    if (hadSelection)
        Roi.paste();
    else
        run("Select None");

    selectImage(duplicateID);
    run("Select None");

    return duplicateID;
}

function assertExpectedImageOrStop(imageID, stageName) {
    selectImage(imageID);
    getDimensions(sw, sh, sc, sz, st);
    if (sw != imgWidth || sh != imgHeight || sc*sz*st != 1)
        exit("Image/window mismatch at " + stageName +
             ". ID=" + imageID + " title=[" + getTitle() + "]" +
             ", actual=" + sw + "x" + sh + " C/Z/T=" + sc + "/" + sz + "/" + st +
             ", expected=" + imgWidth + "x" + imgHeight + " C/Z/T=1/1/1.");
}

function saveTiffOrStop(imageID, path, description) {
    // Always select the exact image ID; do not trust the active Sholl plot window.
    assertExpectedImageOrStop(imageID, "saving " + description + "; " + path);
    saveAs("Tiff", path);
    if (File.exists(path) == 0)
        exit("Saving failed for " + description + ": " + path);
}

function makeWeightedGrayscale(claheID, directionalID, claheWeight, resultTitle) {
    if (claheWeight < 0 || claheWeight > 1)
        exit("CLAHE weight outside [0,1].");
    imageA = duplicateFullImage(claheID, resultTitle + "__temp_clahe");
    selectImage(imageA);
    run("32-bit");
    run("Multiply...", "value=" + claheWeight);

    imageB = duplicateFullImage(directionalID, resultTitle + "__temp_dir");
    selectImage(imageB);
    run("32-bit");
    run("Multiply...", "value=" + (1.0-claheWeight));

    imageCalculator("Add create 32-bit", imageA, imageB);
    fused = getImageID();
    if (fused == imageA || fused == imageB)
        exit("Weighted grayscale blend did not create a distinct image.");
    selectImage(fused);
    getDimensions(fw, fh, fc, fz, ft);
    if (fw != imgWidth || fh != imgHeight || fc*fz*ft != 1 || bitDepth() != 32)
        exit("Weighted grayscale fusion dimensions or bit depth are wrong.");
    rename(resultTitle);
    setMinAndMax(0, 255);
    setOption("ScaleConversions", true);
    run("8-bit");
    if (bitDepth() != 8)
        exit("Weighted grayscale conversion to 8-bit failed.");
    setVoxelSize(PIXEL_WIDTH_UM, PIXEL_HEIGHT_UM, 1.0, ANALYSIS_UNIT);
    selectImage(imageA);
    close();
    selectImage(imageB);
    close();
    selectImage(fused);
    return fused;
}

function makeAutoThresholdMask(grayscaleID, maskTitle, thresholdMethod) {
    binaryID = duplicateFullImage(grayscaleID, maskTitle);
    selectImage(binaryID);
    thresholdOptions = "method=" + thresholdMethod + " white setthreshold";
    if (THRESHOLD_IGNORE_BLACK)
        thresholdOptions = thresholdOptions + " ignore_black";
    if (THRESHOLD_IGNORE_WHITE)
        thresholdOptions = thresholdOptions + " ignore_white";
    run("Auto Threshold", thresholdOptions);
    getThreshold(thresholdMin, thresholdMax);
    if (thresholdMin < 0)
        exit("Auto Threshold failed for " + thresholdMethod);
    setOption("BlackBackground", true);
    run("Convert to Mask");
    selectImage(binaryID);
    rename(maskTitle);
    binaryQcOrStop(maskTitle);
    return newArray(binaryID, thresholdMin);
}

function runDirectionalFilter(
    sourceImageID,
    outputTitle,
    type,
    operation,
    lineLength,
    directions
) {
    selectImage(sourceImageID);

    defaultTitle =
        getTitle() + "-directional";

    imageCountBefore =
        nImages;

    run(
        "Directional Filtering",
        "type=" + type +
        " operation=" + operation +
        " line=" + lineLength +
        " direction=" + directions
    );

    if (nImages != imageCountBefore + 1)
        exit(
            "Directional Filtering did not create " +
            "the expected output."
        );

    if (isOpen(defaultTitle))
        selectWindow(defaultTitle);

    outputID =
        getImageID();

    if (outputID == sourceImageID)
        exit(
            "Directional Filtering output " +
            "could not be identified."
        );

    rename(outputTitle);

    return outputID;
}

function binaryQcOrStop(label) {
    isEightBit =
        bitDepth() == 8;

    isBinaryImage =
        is("binary");

    hasWhite = 0;
    hasBlack = 0;

    if (isEightBit) {
        getHistogram(
            qcValues,
            qcCounts,
            256
        );

        hasWhite =
            qcCounts[255] > 0;

        hasBlack =
            qcCounts[0] > 0;
    }

    print(
        "=== Binary QC: " +
        label + " ==="
    );

    print(
        "8-bit image: " +
        passFail(isEightBit)
    );

    print(
        "Binary image: " +
        passFail(isBinaryImage)
    );

    print(
        "White foreground present: " +
        passFail(hasWhite)
    );

    print(
        "Black background present: " +
        passFail(hasBlack)
    );

    print("==============================");

    if (
        isEightBit == 0 ||
        isBinaryImage == 0 ||
        hasWhite == 0 ||
        hasBlack == 0
    )
        exit(
            "Binary QC failed. " +
            "See the Log window."
        );
}

function passFail(value) {
    if (value)
        return "PASS";

    return "FAIL";
}

function countWhitePixels() {
    getHistogram(
        whiteValues,
        whiteCounts,
        256
    );

    return whiteCounts[255];
}

function askSingleObjectQc() {
    passChoice =
        "PASS - Continue to Sholl";

    editChoice =
        "EDIT & RETRY - Remove unwanted regions";

    Dialog.create(
        "Manual QC of single object"
    );

    Dialog.addMessage(
        "Inspect the retained white object.\n\n" +
        "PASS: continue to Sholl analysis.\n" +
        "EDIT & RETRY: remove unwanted regions " +
        "and extract again.\n\n" +
        "Click Cancel to stop."
    );

    Dialog.addChoice(
        "QC decision:",
        newArray(
            passChoice,
            editChoice
        ),
        passChoice
    );

    Dialog.show();

    return Dialog.getChoice();
}

function manuallyEditSingleObject(
    imageID,
    imageTitle
) {
    selectImage(imageID);
    run("Select None");

    setForegroundColor(0, 0, 0);
    setBackgroundColor(0, 0, 0);
    setTool("freehand");

    pixelsBeforeManualEdit =
        countWhitePixels();

    waitForUser(
        "Edit SingleObject, then click OK",
        "Edit this image:\n" +
        imageTitle + "\n\n" +
        "Remove every unwanted region before clicking OK.\n\n" +
        "Use Freehand/Polygon ROI + Delete repeatedly,\n" +
        "or paint with a black Pencil/Brush.\n\n" +
        "Do not add or reconnect white foreground.\n" +
        "Click OK when editing is finished."
    );

    selectImage(imageID);
    run("Select None");

    pixelsAfterManualEdit =
        countWhitePixels();

    if (
        pixelsAfterManualEdit >
        pixelsBeforeManualEdit
    )
        exit(
            "White foreground pixels were added.\n" +
            "Only black deletion is allowed."
        );

    removedManually =
        pixelsBeforeManualEdit -
        pixelsAfterManualEdit;

    print(
        "Manual edit removed " +
        removedManually +
        " white pixels."
    );

    return removedManually;
}

function retainSomaConnectedComponent(
    imageID,
    centerX,
    centerY,
    preferredSeedX,
    preferredSeedY,
    w,
    h,
    radius
) {
    selectImage(imageID);
    run("Select None");

    seedIsUsable = 0;

    if (
        preferredSeedX >= 0 &&
        preferredSeedX < w &&
        preferredSeedY >= 0 &&
        preferredSeedY < h
    ) {
        dxSeed =
            preferredSeedX -
            round(centerX);

        dySeed =
            preferredSeedY -
            round(centerY);

        seedDistanceSquared =
            dxSeed * dxSeed +
            dySeed * dySeed;

        if (
            seedDistanceSquared <=
            radius * radius &&
            getPixel(
                preferredSeedX,
                preferredSeedY
            ) == 255
        )
            seedIsUsable = 1;
    }

    if (seedIsUsable) {
        cleanupSeedX =
            preferredSeedX;

        cleanupSeedY =
            preferredSeedY;
    } else {
        nearestSeed =
            findNearestWhitePixel(
                centerX,
                centerY,
                w,
                h,
                radius
            );

        if (nearestSeed[3] == 0)
            exit(
                "No white pixel was found within " +
                radius +
                " px of the soma center."
            );

        cleanupSeedX =
            nearestSeed[0];

        cleanupSeedY =
            nearestSeed[1];
    }

    retainedPixels =
        retainWhiteComponent(
            cleanupSeedX,
            cleanupSeedY
        );

    print(
        "Debris cleanup seed: " +
        cleanupSeedX + ", " +
        cleanupSeedY
    );

    print(
        "Soma-connected component: " +
        retainedPixels +
        " white pixels retained"
    );

    return newArray(
        cleanupSeedX,
        cleanupSeedY,
        retainedPixels
    );
}

function findNearestWhitePixel(
    centerX,
    centerY,
    w,
    h,
    radius
) {
    centerPixelX =
        round(centerX);

    centerPixelY =
        round(centerY);

    nearestX = -1;
    nearestY = -1;
    nearestDistanceSquared = 1e30;
    found = 0;

    for (
        nearestDY = -radius;
        nearestDY <= radius;
        nearestDY++
    ) {
        nearestSearchY =
            centerPixelY +
            nearestDY;

        for (
            nearestDX = -radius;
            nearestDX <= radius;
            nearestDX++
        ) {
            nearestSearchX =
                centerPixelX +
                nearestDX;

            currentDistanceSquared =
                nearestDX * nearestDX +
                nearestDY * nearestDY;

            insideImage =
                nearestSearchX >= 0 &&
                nearestSearchX < w &&
                nearestSearchY >= 0 &&
                nearestSearchY < h;

            insideCircle =
                currentDistanceSquared <=
                radius * radius;

            if (
                insideImage &&
                insideCircle &&
                getPixel(
                    nearestSearchX,
                    nearestSearchY
                ) == 255 &&
                currentDistanceSquared <
                nearestDistanceSquared
            ) {
                nearestX =
                    nearestSearchX;

                nearestY =
                    nearestSearchY;

                nearestDistanceSquared =
                    currentDistanceSquared;

                found = 1;
            }
        }
    }

    return newArray(
        nearestX,
        nearestY,
        nearestDistanceSquared,
        found
    );
}

function findLargestComponentNear(
    imageID,
    centerX,
    centerY,
    w,
    h,
    radius,
    scratchTitle
) {
    scratchID =
        duplicateFullImage(
            imageID,
            scratchTitle
        );

    centerPixelX =
        round(centerX);

    centerPixelY =
        round(centerY);

    bestX = -1;
    bestY = -1;
    bestPixels = 0;
    bestDistanceSquared = 1e30;
    found = 0;

    setForegroundColor(0, 0, 0);

    for (
        searchDY = -radius;
        searchDY <= radius;
        searchDY++
    ) {
        searchY =
            centerPixelY +
            searchDY;

        for (
            searchDX = -radius;
            searchDX <= radius;
            searchDX++
        ) {
            searchX =
                centerPixelX +
                searchDX;

            distanceSquared =
                searchDX * searchDX +
                searchDY * searchDY;

            insideImage =
                searchX >= 0 &&
                searchX < w &&
                searchY >= 0 &&
                searchY < h;

            insideCircle =
                distanceSquared <=
                radius * radius;

            if (
                insideImage &&
                insideCircle &&
                getPixel(
                    searchX,
                    searchY
                ) == 255
            ) {
                doWand(
                    searchX,
                    searchY,
                    0,
                    "8-connected"
                );

                if (selectionType != -1) {
                    getHistogram(
                        componentValues,
                        componentCounts,
                        256
                    );

                    componentPixels =
                        componentCounts[255];

                    nearest =
                        nearestPointInSelection(
                            centerX,
                            centerY,
                            w,
                            h,
                            radius
                        );

                    betterSize =
                        componentPixels >
                        bestPixels;

                    sameSizeCloser =
                        componentPixels ==
                        bestPixels &&
                        nearest[2] <
                        bestDistanceSquared;

                    if (
                        nearest[0] >= 0 &&
                        (
                            betterSize ||
                            sameSizeCloser
                        )
                    ) {
                        bestX = nearest[0];
                        bestY = nearest[1];
                        bestPixels =
                            componentPixels;

                        bestDistanceSquared =
                            nearest[2];

                        found = 1;
                    }

                    run("Fill");
                    run("Select None");
                }
            }
        }
    }

    selectImage(scratchID);
    close();

    selectImage(imageID);
    run("Select None");

    return newArray(
        bestX,
        bestY,
        bestPixels,
        found
    );
}

function nearestPointInSelection(
    centerX,
    centerY,
    w,
    h,
    radius
) {
    centerPixelX =
        round(centerX);

    centerPixelY =
        round(centerY);

    nearestX = -1;
    nearestY = -1;
    nearestDistanceSquared = 1e30;

    for (
        nearDY = -radius;
        nearDY <= radius;
        nearDY++
    ) {
        nearY =
            centerPixelY +
            nearDY;

        for (
            nearDX = -radius;
            nearDX <= radius;
            nearDX++
        ) {
            nearX =
                centerPixelX +
                nearDX;

            nearDistanceSquared =
                nearDX * nearDX +
                nearDY * nearDY;

            if (
                nearX >= 0 &&
                nearX < w &&
                nearY >= 0 &&
                nearY < h &&
                nearDistanceSquared <=
                radius * radius &&
                selectionContains(
                    nearX,
                    nearY
                ) &&
                getPixel(
                    nearX,
                    nearY
                ) == 255 &&
                nearDistanceSquared <
                nearestDistanceSquared
            ) {
                nearestX = nearX;
                nearestY = nearY;

                nearestDistanceSquared =
                    nearDistanceSquared;
            }
        }
    }

    return newArray(
        nearestX,
        nearestY,
        nearestDistanceSquared
    );
}

function retainWhiteComponent(
    seedX,
    seedY
) {
    if (getPixel(seedX, seedY) != 255)
        exit(
            "The Wand seed is not white."
        );

    doWand(
        seedX,
        seedY,
        0,
        "8-connected"
    );

    if (selectionType == -1)
        exit(
            "Wand could not select a component."
        );

    setBackgroundColor(0, 0, 0);
    run("Clear Outside");
    run("Select None");

    if (
        bitDepth() != 8 ||
        is("binary") == 0
    )
        exit(
            "The retained object is not " +
            "an 8-bit binary image."
        );

    getHistogram(
        objectValues,
        objectCounts,
        256
    );

    return objectCounts[255];
}

function getImageGeometry(
    centerX,
    centerY,
    w,
    h,
    pixelWidth,
    pixelHeight
) {
    leftDistance =
        centerX * pixelWidth;

    rightDistance =
        (w - centerX) *
        pixelWidth;

    topDistance =
        centerY * pixelHeight;

    bottomDistance =
        (h - centerY) *
        pixelHeight;

    completeRadius =
        min4(
            leftDistance,
            rightDistance,
            topDistance,
            bottomDistance
        );

    corner1 =
        calibratedDistance(
            centerX,
            centerY,
            0,
            0,
            pixelWidth,
            pixelHeight
        );

    corner2 =
        calibratedDistance(
            centerX,
            centerY,
            w,
            0,
            pixelWidth,
            pixelHeight
        );

    corner3 =
        calibratedDistance(
            centerX,
            centerY,
            0,
            h,
            pixelWidth,
            pixelHeight
        );

    corner4 =
        calibratedDistance(
            centerX,
            centerY,
            w,
            h,
            pixelWidth,
            pixelHeight
        );

    farthestCorner =
        max4(
            corner1,
            corner2,
            corner3,
            corner4
        );

    return newArray(
        completeRadius,
        farthestCorner
    );
}

function getObjectMaxRadius(
    centerX,
    centerY,
    w,
    h,
    pixelWidth,
    pixelHeight
) {
    maximumRadius = 0;

    for (
        objectY = 0;
        objectY < h;
        objectY++
    ) {
        for (
            objectX = 0;
            objectX < w;
            objectX++
        ) {
            if (
                getPixel(
                    objectX,
                    objectY
                ) == 255
            ) {
                radius =
                    calibratedDistance(
                        centerX,
                        centerY,
                        objectX,
                        objectY,
                        pixelWidth,
                        pixelHeight
                    );

                if (
                    radius >
                    maximumRadius
                )
                    maximumRadius =
                        radius;
            }
        }
    }

    return maximumRadius;
}

function calibratedDistance(
    x1,
    y1,
    x2,
    y2,
    pixelWidth,
    pixelHeight
) {
    dx =
        (x2 - x1) *
        pixelWidth;

    dy =
        (y2 - y1) *
        pixelHeight;

    return sqrt(
        dx * dx +
        dy * dy
    );
}

function min4(a, b, c, d) {
    minimum = a;

    if (b < minimum)
        minimum = b;

    if (c < minimum)
        minimum = c;

    if (d < minimum)
        minimum = d;

    return minimum;
}

function max4(a, b, c, d) {
    maximum = a;

    if (b > maximum)
        maximum = b;

    if (c > maximum)
        maximum = c;

    if (d > maximum)
        maximum = d;

    return maximum;
}

function whiteTouchesImageEdge(w, h) {
    for (
        edgeX = 0;
        edgeX < w;
        edgeX++
    ) {
        if (
            getPixel(edgeX, 0) == 255 ||
            getPixel(
                edgeX,
                h - 1
            ) == 255
        )
            return 1;
    }

    for (
        edgeY = 0;
        edgeY < h;
        edgeY++
    ) {
        if (
            getPixel(0, edgeY) == 255 ||
            getPixel(
                w - 1,
                edgeY
            ) == 255
        )
            return 1;
    }

    return 0;
}

function setSntShollRadiusPreferences(
    startRadius,
    stepSize,
    endRadius
) {
    script =
        "var P=Packages.java.util.prefs.Preferences;" +
        "var C=Packages.sc.fiji.snt.plugin." +
        "ShollAnalysisImgCmd.class;" +
        "var N=P.userNodeForPackage(C)." +
        "node(C.getSimpleName());" +
        "N.putDouble('startRadius'," +
        startRadius + ");" +
        "N.putDouble('stepSize'," +
        stepSize + ");" +
        "N.putDouble('endRadius'," +
        endRadius + ");" +
        "'ok';";

    preferenceResult =
        eval("script", script);

    if (preferenceResult != "ok")
        exit(
            "Could not initialize " +
            "the Sholl radius settings."
        );
}

function verifyShollRadiiOrStop(
    csvPath,
    expectedStart,
    expectedStep,
    expectedEnd,
    tolerance
) {
    csvText =
        replace(
            File.openAsString(csvPath),
            "\r",
            ""
        );

    csvLines =
        split(csvText, "\n");

    radiusCount = 0;
    firstRadius = -1;
    previousRadius = -1;
    lastRadius = -1;

    for (
        lineIndex = 1;
        lineIndex < csvLines.length;
        lineIndex++
    ) {
        if (lengthOf(csvLines[lineIndex]) > 0) {
            columns =
                split(
                    csvLines[lineIndex],
                    ","
                );

            if (columns.length >= 2) {
                currentRadius =
                    parseFloat(
                        columns[1]
                    );

                if (radiusCount == 0) {
                    firstRadius =
                        currentRadius;
                } else if (
                    abs(
                        (
                            currentRadius -
                            previousRadius
                        ) -
                        expectedStep
                    ) > tolerance
                ) {
                    exit(
                        "Unexpected Sholl radius step.\n\n" +
                        "Expected: " +
                        expectedStep +
                        " analysis units\n" +
                        "Observed: " +
                        (
                            currentRadius -
                            previousRadius
                        ) +
                        " analysis units"
                    );
                }

                previousRadius =
                    currentRadius;

                lastRadius =
                    currentRadius;

                radiusCount++;
            }
        }
    }

    if (radiusCount == 0)
        exit(
            "No sampled Sholl radii were found."
        );

    if (
        abs(
            firstRadius -
            expectedStart
        ) > tolerance
    )
        exit(
            "Unexpected Sholl start radius: " +
            firstRadius + " analysis units."
        );

    if (
        abs(
            lastRadius -
            expectedEnd
        ) > tolerance
    )
        exit(
            "Unexpected Sholl end radius: " +
            lastRadius + " analysis units."
        );
}

function waitForThreeFiles(
    path1,
    path2,
    path3,
    timeoutMs
) {
    elapsedMs = 0;

    while (elapsedMs < timeoutMs) {
        if (
            File.exists(path1) &&
            File.exists(path2) &&
            File.exists(path3)
        )
            return 1;

        wait(250);

        elapsedMs =
            elapsedMs + 250;
    }

    return 0;
}

function renameFileOrStop(
    oldPath,
    newPath,
    description
) {
    if (File.exists(oldPath) == 0)
        exit(
            description +
            " was not created.\n\n" +
            oldPath
        );

    if (File.exists(newPath))
        exit(
            "Do not overwrite: destination appeared " +
            "during analysis.\n\n" +
            newPath
        );

    renameSucceeded =
        File.rename(
            oldPath,
            newPath
        );

    if (
        renameSucceeded == 0 ||
        File.exists(newPath) == 0
    )
        exit(
            "Could not rename " +
            description + ".\n\n" +
            "From: " +
            oldPath + "\n" +
            "To: " +
            newPath
        );
}

function saveLastCsvRowOrStop(
    oldPath,
    newPath,
    description
) {
    if (File.exists(oldPath) == 0)
        exit(
            description +
            " was not created.\n\n" +
            oldPath
        );

    csvText =
        replace(
            File.openAsString(oldPath),
            "\r",
            ""
        );

    csvLines =
        split(
            csvText,
            "\n"
        );

    csvHeader = "";
    lastRow = "";

    for (
        rowIndex = 0;
        rowIndex < csvLines.length;
        rowIndex++
    ) {
        if (lengthOf(csvLines[rowIndex]) > 0) {
            if (csvHeader == "")
                csvHeader =
                    csvLines[rowIndex];
            else
                lastRow =
                    csvLines[rowIndex];
        }
    }

    if (
        csvHeader == "" ||
        lastRow == ""
    )
        exit(
            "Could not identify the current row in " +
            description + "."
        );

    // ---------------- SUMMARY STEP QA (CSV OUTPUT ONLY) ----------------
    // SNT produces a single-row summary that can report 1.9 for a
    // verified profile sampled at 2,4,...,40 um (20 points / 19 gaps).
    // The detailed profile was verified by verifyShollRadiiOrStop()
    // BEFORE this helper runs. Keep SNT's originally reported value in
    // an extra audit column; change ONLY the summary's Radius step field.
    if (indexOf(csvHeader, ",Radius step,") < 0)
        exit("SNT summary header changed; cannot verify Radius step.");

    // Example: Image.tif,"30.351, 30.351, 1",2.0,40.0,1.9,255:255,...
    // Find the comma just after the quoted center. Unlike split(csv,','),
    // this does not confuse commas INSIDE the quoted XYZ center.
    centerEnd = indexOf(lastRow, "\",");
    if (centerEnd < 0)
        exit("Unexpected SNT summary format: quoted center missing.");
    tailStart = centerEnd + 2;
    tail = substring(lastRow, tailStart);

    comma1 = indexOf(tail, ",");
    if (comma1 < 0)
        exit("SNT summary start-radius field missing.");
    rel2 = indexOf(substring(tail, comma1 + 1), ",");
    if (rel2 < 0)
        exit("SNT summary end-radius field missing.");
    comma2 = comma1 + 1 + rel2;
    rel3 = indexOf(substring(tail, comma2 + 1), ",");
    if (rel3 < 0)
        exit("SNT summary radius-step field missing.");
    comma3 = comma2 + 1 + rel3;

    sntStart = parseFloat(substring(tail, 0, comma1));
    sntEnd = parseFloat(substring(tail, comma1 + 1, comma2));
    sntReportedStep = substring(tail, comma2 + 1, comma3);
    if (abs(sntStart - SHOLL_START_UM) > 0.000001 ||
        abs(sntEnd - SHOLL_END_UM) > 0.000001)
        exit("SNT summary start/end does not match verified Sholl settings.");

    verifiedStepStr = d2s(SHOLL_STEP_UM, 6);
    correctedRow = substring(lastRow, 0, tailStart) +
        substring(tail, 0, comma2 + 1) +
        verifiedStepStr + substring(tail, comma3);

    File.saveString(
        csvHeader + ",SNT_Original_Radius_Step,Step_Verification\n" +
        correctedRow + "," + sntReportedStep + ",Detailed_Profile_PASSED\n",
        newPath
    );

    if (File.exists(newPath) == 0)
        exit(
            "Saving failed for " +
            description + "."
        );

    File.delete(oldPath);
}
