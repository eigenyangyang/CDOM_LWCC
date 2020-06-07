Python codes for processing data measured by LWCC (Liquid Waveguide 
Capillary Cells) system. 
Script name: lwcc_postprocess.py

------------
Preparations
------------
Before using the scripts, please:
1 install:
- Python 3.7x or higher version    
- Essential Python packages: numpy, scipy, pandas, matplotlib.

2 correct labelling in .raw files according to the protocol.
Sample labels before the separator (e.g. '#')) need to contain certain keywords:
underway sample - keyword "UW" (case insensitive),
NaCl solution - keyword "NaCl", "PSU" or "Sal" (case insensitive),
and CTD sample labels do not contain any keywords mentioned in the other two catalogs.
 
3 modify the path of the working directory according to your situation in the main program code block (Lines below "if name == 'main':"). 

4 prepare the following files (see examples and data processing procedure below) and put them in the same directory as the one containing .raw files.
- config_PostProc.txt
- rmSpectraIndex.txt
- label_underway.txt and/or label_station.txt
- lwcc_labels_modify.txt (no need if labelling was done as mentioned above.)
- labels_latlon_datetime_lwcc.txt (no need if data are not meant to be uploaded to Pangaea)

------------
Note
------------
The outputs of each data processing step described below are based on the situation where both underway and CTD station samples are available, and that data are meant to be uploaded to Pangaea. If this is not your case, some "outputs" described below will not exist.
For example, if samples are only from CTD stations and the data are not to be uploaded to Pangaea for the moment, after running Step 4.2, the outputs will be
1) a tab delimited absorption data files in the working directory named as
"lwcc_merged_median_ag_station_salcorr_smooth_null_final_expfit.txt" (corrected and fitted absorption spectra by CTD samples) and 
2) a folder named as "ag_station_salcorr_smooth_null_final_expfit" containing plots of corrected and exponentially fitted absorption spectra of CTD samples.

----------------
Data Processing
----------------

-------------------------------------------------------------------
Step 1: Calculate CDOM absorption coefficients from .raw files.
-------------------------------------------------------------------
a) Change the working directory as the one storing .raw data files 
by modifying the variable "wd" (Line ~963).
b) Open the script, comment out Step 2-4 with "#", and then save file.
c) In the terminal, run the script by:
############################################
python3 lwcc_postprocess.py
############################################

Step 1 calculates sample CDOM absorption coefficients by correcting signals from stray light, dark current and reference (MQ). For each measurement, absorption is calculated once, with dark and reference before the sample.

Output: Level 0 absorption coefficients data (.l0a files).

-------------------------------------------------------------------
Step 2: Plot each sample's absorption data in a .l0a file.
-------------------------------------------------------------------
a) Open the script and comment out Step1, 3-4.
b) Run the script.

Output:
Folders named as the base name of the .l0a file in the working direcoty 
containing plots from individual sample measurements. The legend of each 
absorption spectrum in every plot gives the row index of each data point 
when this .l0a data file is imported into a pandas dataframe via "pd.read_csv". 
These row indices can be used as input of Step 3 to remove suspicious data.

--------------------------------------------------------------------------
Step3.1: Manuelly remove suspicious measurements in .l0a data indexed by 
"row_index" (indicated in the plot legend, see Step 2) and take the median 
values and standard deviation of repeatedly measured data.

Step3.2: Extract sample labels from all '*_ag_median.txt' files, merge all '*_ag_median.txt' and '*_ag_sd.txt' files, separate data from CTD staions, underway and NaCl solution, and merge data from each catalog, respectively.
--------------------------------------------------------------------------
a) Open the script and comment out Step1-2, 4.
b) Prepare your own "rmSpectraIndex.txt" file (tab delimited) based on the plots generated from Step2 in the working directory.
c) Run the script.


Output of Step3.1:
1) Tab delimited files in the direcoty "Median_agfiles" named as "*_ag_median.txt" and "*_ag_sd.txt", with "*" standing for the base name of the .l0a file in the working direcoty,respectively.
2) Plots of the spectra in the "*_ag_median.txt" and "*_ag_sd.txt" named as "*_ag_median.png" and "*_ag_sd.png", respectively, in the plotting directory (see Step2). 


Output of Step3.2:
1) Eight tab delimited files in the working directory named as 
    "lwcc_merged_median_ag_salt_uncorrected.txt" (absorption spectra by NaCl solution), 
    "lwcc_merged_sd_ag_salt_uncorrected.txt" (standard deviation of absorption spectra by NaCl solution), 
    "lwcc_merged_median_ag_underway_uncorrected.txt" (uncorrected absorption spectra by underway samples), 
    "lwcc_merged_sd_ag_underway_uncorrected.txt" (standard deviation of absorption spectra by underway samples), 
    "lwcc_merged_median_ag_station_uncorrected.txt" (uncorrected absorption spectra by CTD samples),
    "lwcc_merged_sd_ag_station_uncorrected.txt" (standard deviation of absorption spectra by CTD samples), and
    "lwcc_merged_median_ag_all_uncorrected.txt" (uncorrected absorption spectra by all samples), 
    "lwcc_merged_sd_ag_all_uncorrected.txt" (standard deviation of absorption spectra by all samples), 
    respectively.
2) Two plots named as "lwcc_merged_median_ag_all.png" (all uncorrected absorption spectra)
and "lwcc_merged_sd_ag_all.png" (standard deviation of all absorption spectra),
respectively.
3) A file named as "lwcc_merged_labels.txt" containing all the data lables.
This can be used to construct the files "label_underway.txt", "label_station.txt" and/or "lwcc_labels_modify.txt". 

--------------------------------------------------------------------------
Step 4.1.1: Correct salinity effect and offset in the NIR for LWCC CDOM absorption data (Lefering et al., 2017; Röttgers (2020), in prep). When a spectrum is significantly positive in the NIR (e.g. > 0.05 m^-1 in this case), it is not considered as a valid measurement and removed.

Step 4.1.2: Prepare data suitable for uploading to Pangaea (whether this is needed is specified by the parameter "TFpangaea" in the configuration file "config_PostProc.txt").

Step 4.2: Fit CDOM absorption data with an exponential decay function following Dall’Olmo et al. (2017), i.e. fitting region 420-490 nm, reference wavelength 440 nm.
--------------------------------------------------------------------------
a) Open the script and comment out Step1-3.
b) Prepare your own salinity data for CTD stations and underway data in the working directory, respectively. If they are downloaded from Pangaea, they should be named as "*phys_oce.tab" and "*surf_oce.tab", respectively; if not, rename the files as "*phys_oce_extracted.txt" and "*surf_oce_extracted.txt", respectively.
c) If data are to be uploaded to Pangaea (Step 4.1.2), prepare your own "labels_latlon_datetime_lwcc.txt" in the working directory.; if not, skip this step.
d) Run the script.


Output of 4.1.1:
1) Two tab delimited absorption data files in the working directory named as
"lwcc_merged_median_ag_station_salcorr_smooth_null_final.txt" (corrected absorption spectra by CTD samples) and
"lwcc_merged_median_ag_underway_salcorr_smooth_null_final.txt" (corrected absorption spectra by underway samples),
respectively.
2) Two plots names as
"lwcc_merged_median_ag_station_salcorr_smooth_null_final.png" (corrected  absorption spectra by CTD samples) and
"lwcc_merged_median_ag_underway_salcorr_smooth_null_final.png" (corrected absorption spectra by underway samples),
respectively.


Output of 4.1.2:
1) Six tab delimited absorption data files in the working directory named as
"lwcc_merged_median_ag_station_salcorr_smooth_null_final_pangaea.txt" (corrected absorption spectra by CTD samples),
"lwcc_merged_sd_ag_station_uncorrected_pangaea.txt" (standard deviation of absorption spectra by CTD samples),
"lwcc_merged_median_ag_station_uncorrected_pangaea.txt" (uncorrected absorption spectra by CTD samples),
"lwcc_merged_median_ag_underway_salcorr_smooth_null_final_pangaea.txt" (corrected absorption spectra by underway samples),
"lwcc_merged_median_ag_underway_uncorrected_pangaea.txt" (uncorrected absorption spectra by underway samples), and
"lwcc_merged_sd_ag_underway_uncorrected_pangaea.txt" (standard deviation of absorption spectra by underway samples),
respectively.
2) Two folders named as 
"ag_station_salcorr_smooth_null" and "ag_underway_salcorr_smooth_null", respectively, containing plots of corrected absorption spectra of CTD and underway samples.


Output of 4.2:
1) Four tab delimited absorption data files in the working directory named as
"lwcc_merged_median_ag_station_salcorr_smooth_null_final_expfit.txt" (corrected and fitted absorption spectra by CTD samples),
"lwcc_merged_median_ag_station_salcorr_smooth_null_final_pangaea_expfit.txt" (corrected and fitted absorption spectra by CTD samples),
"lwcc_merged_median_ag_underway_salcorr_smooth_null_final_expfit.txt" (corrected and fitted absorption spectra by underway samples), and
"lwcc_merged_median_ag_underway_salcorr_smooth_null_final_pangaea_expfit.txt" (corrected and fitted absorption spectra by underway samples),
respectively.
2) Four folders named as 
"ag_station_salcorr_smooth_null_final_expfit", "ag_station_salcorr_smooth_null_final_pangaea_expfit", "ag_underway_salcorr_smooth_null_final_expfit", and "ag_underway_salcorr_smooth_null_final_pangaea_expfit", 
respectively, containing plots of corrected and exponentially fitted absorption spectra of CTD and underway samples.


----------
References
----------
1) Lefering et al. (2017). Uncertainty budgets for liquid waveguide CDOM 
absorption measurements. Applied Optics, 56(22), 6357-6366.
2) Röttgers (2020). In prep.
3) Dall’Olmo et al. (2017). Determination of the absorption coefficient of 
chromophoric dissolved organic matter from underway spectrophotometry. Optics 
Express, 25(24), A1079-A1095.
