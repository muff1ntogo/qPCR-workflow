// ============================================================
// GFP Colony Bulk Characterisation Macro
// Fiji/ImageJ — ImageJ Macro Language
// ------------------------------------------------------------
// Measures IntDen, Area, and Mean GFP signal for each
// auto-detected colony across all images in a folder.
// GFP = Channel 2 (green channel)
// ============================================================

// ---------- USER PARAMETERS — edit these ----------
var gaussianSigma   = 2;        // blur radius (px) before thresholding
var thresholdMethod = "Otsu";   // Otsu | RenyiEntropy | Triangle | MaxEntropy
var minColonyArea   = 500;      // px² — ignore ROIs smaller than this
var maxColonyArea   = 999999;   // px² — ignore ROIs larger than this (e.g. debris)
var excludeEdges    = true;     // ignore colonies touching image border
// --------------------------------------------------

// Pick input folder and output file
inputDir   = getDirectory("Select folder containing GFP images");
outputFile = inputDir + "GFP_Colony_Results.csv";

// Initialise results file with header
f = File.open(outputFile);
print(f, "Image,Colony_ID,Area_px2,Mean_GFP,IntDen,RawIntDen");
File.close(f);

// Suppress pop-up dialogs during batch run
setBatchMode(true);
roiManager("reset");

fileList = getFileList(inputDir);

for (i = 0; i < fileList.length; i++) {
    fileName = fileList[i];

    // Process only .tif / .tiff / .czi files
    if (!endsWith(toLowerCase(fileName), ".tif")  &&
        !endsWith(toLowerCase(fileName), ".tiff") &&
        !endsWith(toLowerCase(fileName), ".czi"))  continue;

    open(inputDir + fileName);
    originalID = getImageID();
    imageName  = getTitle();

    // --- 1. Prepare working copy ---
    run("Duplicate...", "title=working");
    workID = getImageID();

    // Split channels and extract GFP (channel 2 / green)
    Stack.getDimensions(w, h, channels, slices, frames);
    if (channels > 1) {
        run("Split Channels");
        selectWindow("C2-working");
        workID = getImageID();
        for (c = 1; c <= channels; c++) {
            if (c != 2 && isOpen("C" + c + "-working"))
                close("C" + c + "-working");
        }
    }

    // Convert to 32-bit for reliable thresholding
    run("32-bit");

    // --- 2. Smooth to reduce noise ---
    run("Gaussian Blur...", "sigma=" + gaussianSigma);

    // --- 3. Auto-threshold → binary mask ---
    setAutoThreshold(thresholdMethod + " dark");
    run("Convert to Mask");

    // Morphological clean-up: fill holes, smooth colony edges
    run("Fill Holes");
    run("Open");   // erode then dilate — removes small speckles

    // --- 4. Detect colonies with Analyze Particles ---
    roiManager("reset");
    if (excludeEdges) {
        edgeFlag = " exclude";
    } else {
        edgeFlag = "";
    }

    run("Analyze Particles...",
        "size=" + minColonyArea + "-" + maxColonyArea +
        " circularity=0.00-1.00" + edgeFlag +
        " add");

    nColonies = roiManager("count");

    // --- 5. Measure IntDen on the ORIGINAL image ---
    selectImage(originalID);

    // Set measurements: area, mean grey, integrated density
    run("Set Measurements...",
        "area mean integrated redirect=None decimal=3");

    for (r = 0; r < nColonies; r++) {
        roiManager("select", r);
        run("Measure");

        // Pull values from the Results table (last row)
        row       = nResults - 1;
        area      = getResult("Area",       row);
        mean      = getResult("Mean",       row);
        intDen    = getResult("IntDen",     row);
        rawIntDen = getResult("RawIntDen",  row);

        // Append to CSV
        File.append(imageName + "," + (r+1) + "," +
                    area + "," + mean + "," +
                    intDen + "," + rawIntDen,
                    outputFile);
    }

    // Tidy up
    close("working");
    selectImage(originalID);
    close();
    roiManager("reset");
    run("Clear Results");

    print("Processed: " + imageName + " — " + nColonies + " colonies found");
}

setBatchMode(false);
print("\nDone. Results saved to:\n" + outputFile);