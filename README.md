## Update 20150217 changes by kmerkens to work with .wav files sampling at rates over 200kHz (no aliasing), particularly for recordings of kogia, dall's porpoise, etc. 
-cat_click_times.m - replaced with cat_click_times_forwav.m
	line 38 "real" times (relative to the baby jesus) calculated from hdr only. 
	line 47 removed lines for calculating "real" times
	line 52 removed lines for calculating "real" times
	line 75 commented out section for plotting - not likely relevant for short wav file samples. 
	line 114 removed save parameters medianValues, meanSpecClicks, iciEncs as not using encounters.
-dLoad_HRsettings
	various parameter thresholds changed for the different example files
-dLoadSTsettings
	various parameter thresholds changed for the different example files
-plotClickEncounters_150310.m
	line 33-36 using clickTimes to calculate ici because they're already dnums relative to the baby jesus, not 		relative to the raw file start like they are in the harp data.
-clickInlinePProc.m
	lines 24-116 commented out code to remove minutes/seconds with too many clicks.
	lines 166-187 aded code to remove echoes immediately following first click for when recording was made in 		captivity.
-dtST_batch.m
	line 29 and 57 added OR statement to deal with formatting of file extension .wav vs. .WAV
-dFind_xwavs.m
	lines 47-48 added code to generate metadata and Graph directories, since it's unlikely that these will ever 		be guided detectors. 
-get_fileset.m
	lines 12-13 added to make.c file names for .wav files, not just .x.wav


## Update 20150205 changes by kmerkens to identify kogia signals on xwav data with 320kHz sampling rate
-Created files: 
	cat_all_mats.m (to combine the .mat files output from cat_click_times.m to get overall summary statistics and 			make plots of the detections from the entire deployment)
	cat_click_times.m (take all the .mat output files from the detector and produce one .mat with all of the 				parameters concatenated for that directory, and to calculate the actual times of the clicks (relative 			to the baby jesus).  Also can generate one set of plots for each encounter. Saves one summary .			mat and one long list of start/end times as .xls
	guidedDetection.m (Added to increase efficiency and accuracy by only running detector over xwav files 			spanned by a previously defined "detection", requires .xls input file, with start/end times of 				encounters formatted as numbers. Called by dFind_xwavs.
	plotClickEncounters_150310.m (Generates plots of clicks according to encounter start/end times, as long as 			the encounter is contained within one .xwav. (so, it's mostly useless))
	plotClickEncounters_posthoc_150310.m (Generates a set of plots for each encounter, even if they 		span multiple xwavs. Called by cat_clicks_mat.m for plotting after the detector has been run.

-de_detector.m
	line 58 added parameter "GuideDetector" to indicate using a pre-definted set of periods ffor checking or not
	line 62 - added input variable "guideDetector" to the dFind_xwavs command, and returned variables 				"encounterTimes" and "GraphDir" to be passed to later children
	line 86 added input variables "encounterTimes,guideDetector,GraphDir" to dHighres_click_batch
-dLoad_HRsettings.m
	various parameters adjusted for different data type
-dLoad_STsettings.m
	various parameters adjusted for different data type
-clickInlinePProc.m
	line 1 added input values "encoutnerTimes, guideDetector,hdr"
	lines 19-69 added code to remove clicks from minutes with more than 75 clicks per minute
	lines 70-116 added code to remove clicks from seconds with more than 20 clicks per second
	lines 161-187 added code to remove clicks from outside encounter times, as identified by guided detector xls 		input sheet. 
-clickParameters.m
	lines 29-30 added parameters zerosvec and adjusted paramter yFilt to be used in lines 57-58
	line 57-58 added clarification to only pull out section of bandwidth allowed through BPF
	line 113 added code to plot the spectrum for troubleshooting (commented out)
	line 150 removed sections for identifying kogia spp in 200 khz data - not applicable. 
-dt)Highres_click_batch.m
	provided input parameters "encounterTimes, guideDetector, GraphDir"
	line 55 added "encounterTimes, guideDetector,hdr" to input parameters for clickInlinePProc.m
	line 82 added call to plotClickEncounters_150310 if guideDetector ==1 to make plots
-dFind_xwavs.m
	provded input parameter guideDetector
	line 10 added additional ffile name to be ignored ("other") in case extra non-xwav files are present (and 		should be ignored) in the main director.
	linee 20-22 added call to guidedDetection if selected


## Update 20150123 changes by kmerkens to identify Kogia spp signals on data from 200kHz sampling rate.
-de_detector.m
	line 25 - added tic
	lines 27-30 added my path to tf file
	line 38  to my basedir
	line 45 to Tinian
	line 54 uncommented indisk…
	lines 73, 74 added “done low res” time elapsed toc
	line 85 added toc
-dLoad_HRsettings.m
	line 45 added 3 parameters for kogia click id, to be used in clickParameters.m
		parametersHR.localminbottom = 63; %Frequency above which will be checked for local min                                         parametersHR.localmintop = 90; %Frequency below which will be checked for local min.                                              parametersHR.localminThr = 68; %Frequency above which a local minimum must exist in order 			to be thrown out.  If the min is below, the click will be kept            
	Multiple parameters changed for Ksp. see comments in file
-dLoad_STsettings.m    
	Multiple parameters changed for Ksp. see comments in file            
- clickParameters
	Throughout - if clicks did not meet a critera, they were eliminated and the code moved on to the next click to 		save time (kfrasier's orinal calculated all parameters, and then checked them at the end)
	line 50 added validClicks = ones(size(ppSignal));
	added sections to do the ffollowing to retain kogia clicks
		-Use the ratio of the median energy near 75 and 97kHz as a cutoff
		-Use the ratio of the median energy near 55 kHz and 70 kHz as a cutoff
		-check for local minimum between 60 and 90 kHz
-dt_highres_click batch
	line 49 commented out message
	line 50 added toc
-dHR_expand_region
	lines 9-14 added check for curve_fitting_toolbox, and use fastsmooth



## Update 11/21/2014
Detector has been updated to accept .wav files in addition to x.wavs.
If you run the detector on directories of wav files, it will look for file start time information in the file name.
 
Edit the regular expression in the load settings scripts:

parametersHR.DateRE = '_(\d*)_(\d*)';

to match your filename date format. 

The result should be a string of numbers in the following order:
yyyymmddHHMMSS

The implementation is a little bit wonky, so contact me if you have problems.

The detector will determine what file type you're using, but to be on the safe side, I'd suggest analyzing wav and xwav files separately.


# Density Estimation Detector

Set of scripts to detect odontocete clicks above a threshold recieved level in .x.wav data. This two-pass detector is based on and incorporates previous work by Marie A. Roch, Simone Baumann-Pickering, and Sean M. Wiggins. 
Detector by Kait E. Frasier.



## How To Use This Code:

1. Place /DE_Detector/ and all subdirectories in Matlab path

2. Edit low resolution detector parameters 
		dLoad_STSettings.m

3. Edit high resolution detector parameters 
		dLoad_HRSettings.m

4. Edit the following lines at the top of *de_detector.m* to reflect your file locations, and which routines to run*

```matlab
% Set transfer function location
tfFullFile = 'E:\Code\TF_files\585_091116_invSensit_MC.tf';
 
% Location of files to be analyzed
baseDir = 'H:\';
 
% Name of the deployment. This should be the first few characters in the
% names of the xwav files you want to look at.
depl = 'GofMX';
 
% Set flags indicating which routines to run. 
lowResDet = 1; %run short time detector.
highResDet = 1; %run high res detector

```

* NOTE: The high res detector relies on output from the low res step, but once you've run the low res, you don't need to re-run it, unless you change your parameters.


## Things to know

 - This detector is set up to run without manual detection inputs. Adding the option to choose times based on a manual pass is trivial, but not currently implemented.

- .wav file inputs are not currently supported.

- Currently optimized for dolphins in the Gulf of Mexico. Use cautiously, especially if adapting to other odontocetes. 

- In noisy areas, a more restrictive postprocessing step is possible. Email for more info.


## Outputs

/metadata/<disk name> 
	contains all of the file types below:
- .c  - 
Low res detector output. This is a text file with flagged start and end times listed as 2 columns. Times are in seconds relative to .x.wav start time.

- .ctg  -   
High res detector output. This is a text file with detection start and end times listed as 2 columns. Times are in seconds relative to .x.wav start time.

- .ptg  - 
same as .ctg but has been run through a post processing step to remove redundant detections, etc.

- .mat  - 
matlab file containing various parameters describing the detected signals retained by all detection steps.

Note: Filenames match the name of the .xwav they describe.


