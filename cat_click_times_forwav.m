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

%Set sampling frequency, in Hz
% fs = 375000; %V Janik
% fs = 500000; %D Mann, E Jacobson
%fs = 480000; %T Yack Dalls
fs = 384000; %CARB

%inDir = 'E:\metadata\bigDL'; % the path to your directory of detector outputs goes here
%inDir = 'D:\metadata\Hawaii18K_disk04';
% inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\VJanik_Ksima_Wild\metadata\kogia';
% inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\DMann_Ksima_captive\metadata\kogia';
%inDir = 'C:\Users\Karlina.Merkens\Documents\Kogia\320_detectctor_dir\metadata\320_Detector_Test';
inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\metadata\kogia';

matList = dir(fullfile(inDir,'kogia*.mat')); % Add wildcard to match the files you want to process.

clickDnum = [];
durClickcon = [];
bw3dbcon = [];
bw10dbcon = [];
nDurcon = [];
ndur95con = [];
ndur95Tailscon = [];
peakFrcon = [];
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
        'dur95','dur95Tails','bw3db','bw10db',...
        'nDur', 'peakFr','ppSignal','specClickTf','specNoiseTf','yFilt','f')
    if ~isempty(clickTimes)
    % determine true click times
        clickDnumTemp = (clickTimes./sec2dnum) + hdr.start.dnum;
        clickDnum = [clickDnum;clickDnumTemp]; %save to one vector
        durClickcon = [durClickcon;durClick];
        bw3dbcon = [bw3dbcon;bw3db];
        bw10dbcon = [bw10dbcon;bw10db];
        nDurcon = [nDurcon; nDur];
        ndur95con = [ndur95con; dur95];
        ndur95Tailscon = [ndur95Tailscon; dur95Tails];
        peakFrcon = [peakFrcon; peakFr];
        ppSignalcon = [ppSignalcon; ppSignal];
        specClickTfcon = [specClickTfcon; specClickTf];
        specNoiseTfcon = [specNoiseTfcon; specNoiseTf];
        yFiltcon = [yFiltcon; yFilt];
        % write label file:
        clickTimeRel = zeros(size(clickDnumTemp));
        % generate label file by replacing .mat extension with .lab for
        % wavesurfer:
        outFileName = strrep(matList(i1).name,'.mat','.lab');
        % open file for writing
        fidOut = fopen(fullfile(inDir,outFileName),'w+');
        fclose(fidOut);
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


%%%Section added to do post-processing where all the clicks are together,
%%%not speparted by xwav. 

%Get detectionTimes
% inpath = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\VJanik_Ksima_Wild\kogia';
% infile = 'VJanik_Ksima_Wild_log_150521.xls';
% inpath = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\DMann_Ksima_captive\kogia';
% infile = 'DMann_Ksima_captive_log_150626.xls';
inpath = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\kogia';
infile = 'Ksima_guided_detector_160601.xls';

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
specClickTf = specClickTfcon;
specNoiseTf = specNoiseTfcon;
peakFr = peakFrcon;
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
    durClickconP, ndur95conP, ndur95TailsconP, bw3dbconP, bw10dbconP, nDurconP, peakFrconP, ppSignalconP,...
    specClickTfconP,specNoiseTfconP, yFiltconP] = plotClickEncounters_posthoc_150310(encounterTimes,...
    clickTimes,ppSignal,durClick,ndur95,ndur95Tails,bw3db,bw10db,...
    specClickTf,specNoiseTf,peakFr,nDur,yFilt,hdr,GraphDir,f);

%Change the name on the pruned parameters
clickDnum = clickTimesconP;
durClickcon = durClickconP;
bw3dbcon = bw3dbconP;
bw10dbcon = bw10dbconP;
nDurcon = nDurconP;
ndur95con = ndur95conP;
ndur95Tailscon = ndur95TailsconP;
peakFrcon = peakFrconP;
ppSignalcon = ppSignalconP;
specClickTfcon = specClickTfconP;
specNoiseTfcon = specNoiseTfconP;
yFiltcon = yFiltconP;

%Then save everything
save([inDir,'\',choppedDir{7},'_ClicksOnlyConcat',filedate,'.mat'],...
    'clickDnum','durClickcon','ndur95con','ndur95Tailscon','bw3dbcon',...
    'bw10dbcon','nDurcon', 'peakFrcon','ppSignalcon',...
    'specClickTfcon','specNoiseTfcon','yFiltcon','medianValues',...
    'meanSpecClicks','meanSpecNoises','iciEncs','f')

%Save the pruned clicks that remain after removing any with too small icis




