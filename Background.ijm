// ============================================================
// Background Calibration Macro
// Run this ONCE on your folder of no-GFP control images
// ============================================================

controlDir = getDirectory("Select folder containing control images");
outputFile = controlDir + "background_value.txt";

fileList = getFileList(controlDir);
setBatchMode(true);

totalBackground = 0;
imageCount = 0;

for (i = 0; i < fileList.length; i++) {
    fileName = fileList[i];

    if (!endsWith(toLowerCase(fileName), ".czi")  &&
        !endsWith(toLowerCase(fileName), ".tiff") &&
        !endsWith(toLowerCase(fileName), ".png"))  continue;

    open(controlDir + fileName);
    originalID = getImageID();

    // Extract channel 1 if multichannel
    Stack.getDimensions(w, h, channels, slices, frames);
    if (channels > 1) {
        run("Split Channels");
        selectWindow("C1-" + fileName);
        for (c = 2; c <= channels; c++) {
            if (isOpen("C" + c + "-" + fileName))
                close("C" + c + "-" + fileName);
        }
    }

    run("32-bit");

    // Measure mean intensity of the entire image as background
    run("Set Measurements...", "mean redirect=None decimal=3");
    run("Select All");
    run("Measure");

    bgValue = getResult("Mean", nResults - 1);
    totalBackground += bgValue;
    imageCount++;

    print("Control image: " + fileName + " | Mean: " + bgValue);

    run("Clear Results");
    close();
}

setBatchMode(false);

// Calculate and save averaged background
avgBackground = totalBackground / imageCount;
print("\nAveraged background across " + imageCount + " controls: " + avgBackground);

// Save to file so the main macro can read it
f = File.open(outputFile);
print(f, avgBackground);
File.close(f);

print("Background value saved to: " + outputFile);