# How to use proper velocity inlets for your (non-cylindrical) geometry:

## Step 1:
Execute the script to generate the weights file itself
python <base_dir>/hemelb-dev/hemelb/Tools/setuptool/InletProcessing/CreateWeightsFile.py <name of the profile file in your configuration directory>

NOTE: The script will also look for your XML and GMY files, so those should reside in the locations indicated in the profile file.

## Step 2:
Ensure HemeLB is properly configured. This means:
a. Option HEMELB_USE_VELOCITY_WEIGHTS_FILE is set to ON
b. You are using Ladd velocity inlet conditions.
c. You are using 0th order pressure outlet conditions.
d. Set HEMELB_COMPUTE_ARCHITECTURE to ISBFILEVELOCITYINLET. This is the weighting configuration for Intel Sandy Bridge, modfied with higher weightings for the inlets. Enabling this option will gain you up to 64% more performance compated to the neutral setting.

NOTE:
The weighting mechanism could work with outlets, but the current CreateWeightsFile.py only generates weighting files for the inlets.

## Step 3:
Compile the code with the new settings.

## Step 4:
Double check your XML file, and make sure you put in all the Property Extraction Planes.

## Step 5:
Run HemeLB :).
