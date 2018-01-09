% cat_click_times.m Code to take all the .mat output files from the
% detector and produce one .mat with all of the parameters concatenated for
% that directory, and to calculate the actual times of the clicks (relative
% to the baby jesus).  Also can call the plotting code to generate one set
% of plots for each encounter (not separated by .xwav file), unless the
% bottom section is commented out.
% Saves one summary .mat and one long list of start/end times as .xls (and
% plots, if requested)
% Different from cat_click_times.m - doesn't calculate time relative to raw
% start, uses hdr info. plotting section commented out.



% ****ALSO SET PATHS FOR GUIDED ANALYSIS BELOW, APPX LINES 120

close all
% clear all

%inDir = 'E:\metadata\bigDL'; % the path to your directory of detector outputs goes here
%inDir = 'D:\metadata\Hawaii18K_disk04';
% inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\VJanik_Ksima_Wild\metadata\kogia';
% inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\DMann_Ksima_captive\metadata\kogia';
%inDir = 'C:\Users\Karlina.Merkens\Documents\Kogia\320_detectctor_dir\metadata\320_Detector_Test';
% inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\metadata\kogia';
% inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\TGridley_Ksima_Wild\metadata\kogia';
inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_DASPR_2017\metadata\kogia';


% matList = dir(fullfile(inDir,'kogia*.mat')); % Add wildcard to match the files you want to process.
matList = dir(fullfile(inDir,'1*.mat')); %for DASPRs it's the first few digits of the instrument name. 
%If you don't speficy this, it will attempt to use any/all .mat files, and
%that can be a problem if you already have run this and have some summary
%.mat files in the directory. 

%ID whether this is a daspr 1 = yes, 0 = no. Determines how bouts are
%output.
DASPR = 1;


clickDnum = [];
durClickcon = [];
bw3dbcon = [];
bw10dbcon = [];
bwRMScon = [];
nDurcon = [];
ndur95con = [];
ndur95Tailscon = [];
peakFrcon = [];
centFrcon = [];
QRMScon = [];
Q3dBcon = [];
snrcon = [];
ppSignalcon = [];
specClickTfcon = [];
specNoiseTfcon = [];
yFiltcon = [];

sec2dnum = 60*60*24; % conversion factor to get from seconds to matlab datenum
% iterate over detector-derived mat files in directory
for i1 = 1:length(matList)
    clickTimes = [];
    clickDnumTemp = [];
    % only need to load hdr and click times
    load(fullfile(inDir,matList(i1).name),'hdr','clickTimes', 'durClick', ...
        'dur95','dur95Tails','bw3db','bw10db','bwRMS','QRMS','Q3dB',...
        'nDur', 'peakFr','centFr','snr','ppSignal','specClickTf','specNoiseTf','yFilt','f')
    if ~isempty(clickTimes)
    % determine true click times
        clickDnumTemp = (clickTimes./sec2dnum) + hdr.start.dnum;
        clickDnum = [clickDnum;clickDnumTemp]; %save to one vector
        durClickcon = [durClickcon;durClick];
        bw3dbcon = [bw3dbcon;bw3db];
        bw10dbcon = [bw10dbcon;bw10db];
        bwRMScon = [bwRMScon;bwRMS];
        nDurcon = [nDurcon; nDur];
        ndur95con = [ndur95con; dur95];
        ndur95Tailscon = [ndur95Tailscon; dur95Tails];
        peakFrcon = [peakFrcon; peakFr];
        centFrcon = [centFrcon; centFr];
        QRMScon = [QRMScon; QRMS];
        Q3dBcon = [Q3dBcon; Q3dB];
        snrcon = [snrcon;snr];
        ppSignalcon = [ppSignalcon; ppSignal];
        specClickTfcon = [specClickTfcon; specClickTf];
        specNoiseTfcon = [specNoiseTfcon; specNoiseTf];
        yFiltcon = [yFiltcon; yFilt];
%         %%%This part isn't working right. 
%         % write label file:use clickTimes, which is relative to the start
%         % of the file, which is what we want for wavesurfer. Add a label:
%         numclicks = size(clickTimes,1);
%         clickLabelsNum = [001:1:numclicks]';
%         clickLabelsStr = num2str(clickLabelsNum);
%         %now make into cell arrays and combine
% %         clickTimesCell = num2cell(clickTimes); 
% %         clickTimesCell = cell(clickTimes);
%          clickLabelsCell = cellstr(clickLabelsStr);
% %         clickLabels = {clickTimesCell,clickLabelsCell};
%         %turn it into a table
% %         clickLabelsTab = cell2table(clickLabels);
%         clickLabelsTab = table(clickTimes,clickLabelsCell);
%         % generate label file by replacing .mat extension with .lab for
%         % wavesurfer:
%         outFileName = strrep(matList(i1).name,'.mat','.txt');
%         outFilePath = [inDir,'\',outFileName];
%         %Write file
%         %csvwrite(outFilePath,clickLabels);
%         writetable(clickLabelsTab,outFilePath)
    end
end

choppedDir = strsplit(inDir,'\'); %cut up the file path to get the disk name
%so that you can save the files with identification. 
filedate = datestr(now, 'yymmdd');


%Added to save start/end times as character arrays
clickDnumChar1 = char(datestr(clickDnum(:,1)));
clickDnumChar2 = char(datestr(clickDnum(:,2)));
numclicks = size(clickDnum,1);
clickDnumChar = {};
for nc = 1:numclicks
    clickDnumChar{nc,1} = clickDnumChar1(nc,:);
    clickDnumChar{nc,2} = clickDnumChar2(nc,:);
end
xlswrite([inDir,'\',choppedDir{7},'_ClicksOnlyConcatCHAR',filedate,'.xls'],clickDnumChar)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Add step to go through and make a "short list" of the detections, using
%clicks separated by not more than 3 minutes, providing a start time, end
%time and the total number of clicks. This is used for verifying the click
%bouts in triton to make a log before running a guided detector.
%180105-changed to be only 2 minutes, to work with the daspr files. Also 
%add file name, for ease of identification
startclick = clickDnum(1,1);
threemin = datenum([0,0,0,0,2,0]); %actually two minute now. 
clickitr = 1;
boutitr = 1;
bouts = [];
for ncc = 2:numclicks
   prevclick = clickDnum(ncc-1,2);
   checkclick = clickDnum(ncc,1);
   clickdiff = checkclick-prevclick;
   if clickdiff > threemin || ncc == numclicks
       bouts(boutitr,1) = startclick;
       if ncc == numclicks
           endclick = clickDnum(ncc,2);
       else
           endclick = clickDnum(ncc-1,2);
       end
       bouts(boutitr,2) = endclick;
       if ncc == numclicks
           clickitr = clickitr +1;
       end
       bouts(boutitr,3) = clickitr;
       if ncc < numclicks
           startclick = clickDnum(ncc,2);
       else
           continue
       end   
       clickitr = 1;
       boutitr = boutitr + 1;
   elseif clickdiff < threemin
       clickitr = clickitr + 1;
       continue
   end
end
boutsChar1 = char(datestr(bouts(:,1)));
boutsChar2 = char(datestr(bouts(:,2)));
boutsChar3 = num2str(bouts(:,3));
numbouts = size(bouts,1);
boutsChar = {};
if DASPR == 1
    %Format the start date/time to match the file names, for easier
    %identification of the correct file. (it's too hard to try to match
    %files)
    boutsChar4 = [];
    for bb = 1:numbouts
        datestartvec = datevec(bouts(bb,1));
        datestartdatemin = datestr(datestartvec,'yymmddHHMM');
        boutsChar4 = [boutsChar4;datestartdatemin];
    end
end
    
for nb = 1:numbouts
    boutsChar{nb,1} = boutsChar1(nb,:);
    boutsChar{nb,2} = boutsChar2(nb,:);
    boutsChar{nb,3} = boutsChar3(nb,:);
    if DASPR == 1
        boutsChar{nb,4} = boutsChar4(nb,:);
    end
end
xlswrite([inDir,'\',choppedDir{7},'_BOUTS',filedate,'.xls'],boutsChar)










% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Section added to do post-processing where all the clicks are together,
%%%not speparted by xwav. 

%Get detectionTimes
% inpath = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\VJanik_Ksima_Wild\kogia';
% infile = 'VJanik_Ksima_Wild_log_150521.xls';
% inpath = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\DMann_Ksima_captive\kogia';
% infile = 'DMann_Ksima_captive_log_150626.xls';
% inpath = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\kogia';
% infile = 'Ksima_guided_detector_160601.xls';
% infile = 'Ksima_guided_detector_161013_noBUZZ.xls';
% infile = 'Ksima_guided_detector_160601_BUZZes.xls';
% inpath = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\TGridley_Ksima_Wild\kogia';
% infile = 'TGridley_Ksima_Wild_log_170511.xls';
inpath = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_DASPR_2017\kogia';
infile = 'DASPR_Kspp_Wild_log_180104.xls';

%read the file into 3 matrices-- numeric, text, and raw cell array
[num, txt, raw] = xlsread([inpath '\' infile]);
%error check
[~,y]=size(num);
if y < 2;          %start and end dates not formatted as numbers
    h=errordlg('Please save dates in number format and click ok');
    uiwait(h)
    [num, txt, raw] = xlsread([inpath '\' infile]); %reread file
end  
excelDates = num(:,1:2);                %numeric array contains datenums
%convert excel datenums to matlab datenums (different pivot year)
matlabDates = ones(size(excelDates)).*datenum('30-Dec-1899') ...
    + excelDates;

%Use other code to do the plotting and get the medians.
%rename things
encounterTimes = matlabDates;
clickTimes = clickDnum;
guideDetector = 1;
ppSignal = ppSignalcon;
durClick = durClickcon;
bw3db = bw3dbcon;
bw10db = bw10dbcon;
bwRMS = bwRMScon;
QRMS = QRMScon;
Q3dB = Q3dBcon;
specClickTf = specClickTfcon;
specNoiseTf = specNoiseTfcon;
peakFr = peakFrcon;
centFr = centFrcon;
snr = snrcon;
nDur = nDurcon;
ndur95 = ndur95con;
ndur95Tails = ndur95Tailscon;
yFilt = yFiltcon;
GraphDir = [inDir,'\matlab_graphs'];

% %This one doesn't output the pruned parameters
% [medianValues,meanSpecClicks,meanSpecNoises,iciEncs] = plotClickEncounters_posthoc_150310(encounterTimes,...
%     clickTimes,ppSignal,durClick,bw3db,bw10db,...
%     specClickTf,specNoiseTf,peakFr,nDur,yFilt,hdr,GraphDir,f);

[medianValues,meanSpecClicks,meanSpecNoises,iciEncs,clickTimesconP,...
    durClickconP, ndur95conP, ndur95TailsconP, bw3dbconP, bw10dbconP, ...
    bwRMSconP, QRMSconP, Q3dBconP, nDurconP, peakFrconP, centFrconP, ...
    snrconP, ppSignalconP, specClickTfconP,specNoiseTfconP,...
    yFiltconP] = plotClickEncounters_posthoc_150310(encounterTimes,...
    clickTimes,ppSignal,durClick,ndur95,ndur95Tails,bw3db,bw10db,bwRMS,QRMS,Q3dB,...
    specClickTf,specNoiseTf,peakFr,centFr,snr,nDur,yFilt,hdr,GraphDir,f);

%Change the name on the pruned parameters
clickDnum = clickTimesconP;
durClickcon = durClickconP;
bw3dbcon = bw3dbconP;
bw10dbcon = bw10dbconP;
bwRMScon = bwRMSconP;
QRMScon = QRMSconP;
Q3dBcon = Q3dBconP;
nDurcon = nDurconP;
ndur95con = ndur95conP;
ndur95Tailscon = ndur95TailsconP;
peakFrcon = peakFrconP;
centFrcon = centFrconP;
ppSignalcon = ppSignalconP;
specClickTfcon = specClickTfconP;
specNoiseTfcon = specNoiseTfconP;
yFiltcon = yFiltconP;
snrcon = snrconP;

%Then save everything
save([inDir,'\',choppedDir{7},'_ClicksOnlyConcat',filedate,'.mat'],...
    'clickDnum','durClickcon','ndur95con','ndur95Tailscon','bw3dbcon',...
    'bw10dbcon','bwRMScon','QRMScon','Q3dBcon','nDurcon', 'peakFrcon',...
    'centFrcon','snrcon','ppSignalcon',...
    'specClickTfcon','specNoiseTfcon','yFiltcon','medianValues',...
    'meanSpecClicks','meanSpecNoises','iciEncs','f','hdr')

%Save the pruned clicks that remain after removing any with too small icis




