function de_detector
% This is the starting point for a simplified detector based on Marie Roche's  
% teager energy detector. It includes ideas/code snips from Simone, 
% and calls functions from triton. 
% The goal of this detector is to have predictable performance, 
% for use with model-based density estimation efforts. To accomplish this,
% it uses a simple energy threshold to identify clicks, thereby reducing
% the impact of changing noise conditions on detectability. 

% Known issue: Prop cavitation noise often makes it through detector and
% classifier steps.

% The low and hi-res detection passes still happen, but no teager energy 
% is used.

% All input parameters are contained within two separate scripts:
%   dLoad_STsettings : settings for low res detector
%   dLoad_HRsettings : settings for hi res detector
% See those files for info on settings.

% NOTE: If this is a captive recording, there is a "lockout" period that
% can be set in clickInlinePProc.m, line 161 to remove echoes

% clearvars
close all
fclose all;
warning('off','all') %This useful in the more recent versions of matlab because there are a LOT Of warnings.
tic

% Set transfer function location
%tfFullFile = 'E:\Code\TF_files\604_100614\604_100614_invSensit.tf';
%tfFullFile = 'H:\Cetacean Research Program\HARP\TF_files\695_121203_invSensit.tf';
%tfFullFile = 'C:\Users\Karlina.Merkens\Documents\HARPTFfiles\600_series\692_121119\692_121119_invSensit.tf';
%tfFullFile = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\CARB_TF.tf';

%Note, if you don't have a tranfer function just use:
tfFullFile = [];


% Location of base directory containing directories of files to be analyzed
%baseDir = 'I:\GofMXArraySpRecs\Sb';
%baseDir = 'D:\';
% baseDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\VJanik_Ksima_Wild\';
% baseDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\DMann_Ksima_captive\';
% baseDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\';
% baseDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\TGridley_Ksima_Wild\';
% baseDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_DASBR_2017\';
baseDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\Jenner_Ksima_Wild\';

% Name of the deployment. This should be the first few characters in the 
% directory(ies) you want to look in you want to look at. For now,
% directory hierarchy is expected to be: basedir>depl*>*.x.wav
% TODO: implement recursive directory search for more flexibility.
%depl = 'Hawaii';

depl = 'Kogia';
% depl = 'kogia';
%depl = 'dalls';
%depl = 'harbor';

DASBR = 0; %Change this to be 0 if this is not a daspr file. 
%If it IS, then make this 1, and the code will ignore any clicks with a
%time stamp of 0.06 or less, which seem to appear at the start of all
%files. 

% Set flags indicating which routines to run. 
lowResDet = 0; %run short time detector.
highResDet = 1; %run high res detector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[metaDir] = dBuild_dirs(baseDir);
inDisk = fileparts(baseDir(1:3));

% Build list of (x)wav names in the base directory.
% Right now only wav and xwav files are looked for.
guideDetector = 0; %1 if using xls sheet to guide detection, 0 to run on all files in drive

[detFiles,encounterTimes,GraphDir]= dFind_xwavs(baseDir,depl,guideDetector); 

viewPath = {metaDir, baseDir};
[fullFiles,fullLabels] = get_fileset(baseDir,metaDir,detFiles); % returns a list of files to scan through
% profile on
% profile clear
if ~isempty(detFiles)
    % Short time detector
    if lowResDet == 1
        % load settings
        parametersST = dLoad_STsettings;
        % run detector
        dtST_batch(baseDir,detFiles,parametersST,viewPath);
        display('Done with low-res detector')
        toc
    end
    
    % High res detector
    if highResDet == 1
        % load settings
        parametersHR = dLoad_HRsettings;
        % run detector
        dHighres_click_batch(fullFiles,fullLabels,baseDir,parametersHR,...
            viewPath,tfFullFile,encounterTimes,guideDetector,GraphDir,DASBR)
        display('Done with high-res detector')
    end
end

toc
% profile viewer
% profile off