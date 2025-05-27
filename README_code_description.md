MATLAB Batch Processing of TEMPO Scans
A MATLAB package to generate clean and registered data for each scan of TEMPO Level 1 twilight radiance files for a given time range.    The data files provide radiance measured during twilight hours to capture city lights at TEMPO’s native spatial resolution, ~10 km^2 at the center of the Field of Regard (FOR), for individual granules. In this package the processing is performed on a month-worth of data at a time.  The TEMPO RADT granules are downloaded from NASA Earthdata, https://search.earthdata.nasa.gov/search
# Description
The batch scan processing package performs the following:
1.	Loops over all granules in a directory and finds the first granule of a scan to use for extracting two patterns, a time stamp and a scan number, that are used to find all the granules of a scan.
2.	Processing the granules of each includes the following:
a.	Cleaning the VIS and UV images; this includes correction of any speckling and streaks and background removal. 
b.	Saving the cleaned granule radiance and related data to MATLAB files; VIS and UV are saved separately. 
c.	Gathering ephemeris information to be used in computation of Sun and Moon illumination and viewing geometry
d.	Integration of radiance over spectral ranges to form radiance images
e.	Registration of the clean integrated images to a VIIRS DNB clear-sky mosaic
f.	Saving of integrated radiances on a canvas which is the combination of the granules of a scan
g.	Detection of lightning and city lights
h.	Saving of scan variables to a MATLAB file for each scan
2.	Requirements
- MATLAB R2024b or newer (older versions may work but are untested)
- Parallel processing toolbox
3.	Installation
 Clone the repository below or download the files to your local directory.
https://github.com/Smithsonian/TEMPO-MATLAB-processing
# Usage
In MATLAB command line open and run Batch_Process_Scans.
Notes:  In Batch_Process_Scans we need to specify the directory where the raw TEMPO granules are (p_in) and the output directory where the cleaned granules, processed scans, and saved figures will be saved (p_out).
We need to also specify the start and stop day of the data to process. An example is shown below for the data of the month of October 2024.
p_in  =  ‘H:\TEMPO_DATA_oct2024’
p_out  = ‘H:\TEMPO_DATA_oct2024\Cleaned’
start = 20241001;   % yyyymmdd
stop = 20241031;
Note that in this version of the code the Cleaned directory must exist before running the process.  The subdirectories under Cleaned are created by the process.

 
# Files for the SCANS Batch Processing process
| File	| Description |
| Batch_Process_Scans.m |	Main script to process TEMPO scans.  Works with TEMPO granules that are stored in a directory.  Gets the listing of the directory and loops over all the TEMPO granules in that directory to find the granule number from the file name. It then loops over all the found granule numbers and skip them unless it’s granule 1 or the first granule of a scan.  In this latter case a pattern consisting of a time stamp and a scan number are extracted from the granule filename.  The pattern is of the form:  {yyyymmdd, Sxxx}.  The time stamp is used to verify if the granule is for a pre-selected time range and in that case the function Process_Scan is called to process all the granules of the scan associated with the selected pattern.|
|Configvars_Process_Scan.m|	Configuration file contains:
Clean up flags [vis uv]
VIIRS-DNB registration
Standard fixed grid sampling distance in rad
Compositing with moonlight flag
Parameters for city lights classifier|
| Process_Scan.m	| This function is called by the main script Batch_Process_Scan for each scan.
Runs the configuration file
Gets files listing from input directory
Finds RADT netCDF files for granules belonging to scan and for each granule:
Gets ephemeris data from a granule and interpolates satellite ephemeris to pixel times.
Repairs latitude and longitude to fix telemetry defects.
Processes VIS data and then UV data by performing clean up and background removal (if flags enabled).
Converts radiance from units with photons/s to ones with nW.
Integrates over spectral ranges to form radiance image.
Paints granule radiance onto scan canvas.
for UV, integration over lightning spectral range is performed to form radiance image
Saves cleaned granule radiance and related data; one file for VIS and one for UV.
Once the steps above are completed then more processing is performed for the scan:
Subtracts biases that are due to calibration and background subtraction artifacts using the debias function.
Makes figures,  collects statistics and saves with figures.
Registers images to a VIIRS DNB clear-sky mosaic.  This includes computation of the geometric error metrics relative to VIIRS-DNB.
Extends TEMPO navigation beyond limb.
Computes Lunar illumination and satellite viewing geometry.
Saves all relevant data into a scan file.|
| repair_lat_lon.m |	Function to remove pixels with bad quality flag and to perform linear extrapolation of the latitude and longitude to good quality pixels.|
| clean_up.m	| Function to clean the radiance of a granule.  Applies to VIS and UV.
Finds anomalously oversubtracted quadrants.
Finds persistently hot pixels.
Performs despeckle/destreak  on each frame. 
Assigns quality flag value to bad quadrants.| 
| despeckle.m	| Function to find persistently hot pixels in TEMPO radiance frames| 
| despeckle2.m	| Function to clean radiance frames from specks.| 
| destreak2	| Function to clean radiance frames from streaks.| 

| model_bkgnd4.m| 	Function to compute and remove background from radiance images.
Integrates over spectral bins. 
Finds bright locations (cities) in top and bottom of frame using a  MAD filter.
Fits a background for each spectral bin and removes the background from the signal.  Background and signal are returned for saving in the scan file.| 
| bucket.m	| Function to integrate radiance over light bucket (wavelength range)| 
| output_L1p5_granule_file.m	|  to save cleaned radiances for each granule.  VIS and UV data are saved in separate files. | 
| debias.m	| Function to subtract from cleaned radiances biases that are due to calibration and background subtraction artifacts| 
| Make_ref_viirs_hires.m	| Function to make a higher resolution VIIRS-DNB image.| 
| block_bin.m	| Function for binning imagery when needed.| 
| align_tile.m	| Function used in registration between TEMPO and VIIRS-DNB| 
| refine3a.m	| Function used to refine the registration between TEMPO and VIIRS DNB at subpixel resolution.| 
| fill_nan_holes| 	Function to fill double line holes between tiles where overlap was insufficient| 
| madfilt.m	| Function for Median Absolute Difference (MAD) filtering| 
| fixed_grid.m	| Function to determine the navigation variables for TEMPO fixed grid| 
| extend_fg_coords.m	| Function to extend TEMPO fixed grid coordinates beyond the limb| 
| inv_fixed_grid.m	| Function to get latitude, longitude, and height of a site from fixed grid coordinates.| 
| SiteCoordinates.m	| Function to compute site coordinates of point on earth's surface| 
| Pierce2.m	| Function to locate point where line-of-sight pierces the ellipsoid| 
| show_sun_moon.m	| Function to display Earth, sun, moon and terminators| 
| terminator.m	| Function to determine a terminator curve on a spherical Earth illuminated by s (s is the Sun or Moon)| 
| output_L2_scan_files.m	| Function to output to a file all the variables associated to a processed scan.| 
 
Once scan data has been obtained after running Batch_Process_Scans it is possible to generate a clear sky mosaic using the Make_clear_sky_mosaic.m script and to classify the radiance with a spectral library using the Scan_classifier.m script.  These processes are described here. 
Make_clear_sky_mosaic
1.	Uses the saved Level 2 data to make a mosaic of clear sky for a whole month. 
2.	Compares TEMPO data to VIIRS-DNB.  
The files in the table below are used in Make_clear_sky_mosaic.m
| File	| Description| 
| Make_clear_sky_mosaic.m	| Main script for generating a clear sky mosaic of cleaned TEMPO scans.  The script loops over available scans and 
Skips scans that were not completely cleaned and registered to VIIRS.
Masks out moonlit pixels.
Masks lightning.
Generates the mosaic. 
Checks the registration between TEMPO and VIIRS-DNB.
Performs a radiance comparison between TEMPO and VIIRS-DNB.| 
| ligntning_detection.m	| Function to detect pixels where lightning is present and to set those pixels to -inf to mask them.| 
| inv_fixed_grid.m	| Function to get latitude, longitude, and height of a site from fixed grid coordinates.| 
| Make_ref_viirs_hires.m| 	Function to make a higher resolution VIIRS-DNB image.| 
| block_bin.m	| Function for binning imagery when needed.| 
| align_tile.m	| Function used in registration between TEMPO and VIIRS-DNB.| 
| density_histogram.m	| Function to generate density histograms of TEMPO radiance versus VIIRS-DNB radiance.| 
 

Scan_classifier
Goes over each scan file and for each cleaned granule of that scan performs a classification of the radiance using a spectral library


| File	| Description| 
| Scan_classifier.m	Main|  script for classification of TEMPO city lights using a spectral light library.  
The script loops over available scans, extracts the {timestamp scan number} pattern from the scan file name and calls the classification function Classify_Scan.m. | 

| Classify_Scan.m| 	Function to perform classification of TEMPO city lights to a spectral light library.  
Uses the pattern from the main script to find all the granules belonging to the scan.
Identifies the categories of the spectral library.
Normalizes each library spectrum for unit response = 1 nW/(cm^2 sr) across the appropriate spectral range (for VIS or UV) and converts units from photon to energy.
Loops over the granules belonging to the scan and classifies bright pixels.  | 

| Stepwise_Regression_Classification.m	| Function to perform stepwise regression classification of the TEMPO radiance.| 
| wght_spectral_fit | Function to perform a weighted spectral fit.| 

