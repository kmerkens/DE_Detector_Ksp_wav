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
%fs = 375000; %V Janik
fs = 500000; %D Mann, E Jacobson
%fs = 480000; %T Yack Dalls

%inDir = 'E:\metadata\bigDL'; % the path to your directory of detector outputs goes here
%inDir = 'D:\metadata\Hawaii18K_disk04';
inDir = 'C:\Users\Karlina.Merkens\Documents\Porpoise\OtherRecordings\EJacobson_Harbor_wild\metadata\harbor';
%inDir = 'C:\Users\Karlina.Merkens\Documents\Kogia\320_detectctor_dir\metadata\320_Detector_Test';
matList = dir(fullfile(inDir,'harbor*.mat')); % Add wildcard to match the files you want to process.
clickDnum = [];
durClickcon = [];
nDurcon = [];
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
        'nDur', 'peakFr','ppSignal','specClickTf','specNoiseTf','yFilt','f')
    if ~isempty(clickTimes)
    % determine true click times
        clickDnumTemp = (clickTimes./sec2dnum) + hdr.start.dnum;
        clickDnum = [clickDnum;clickDnumTemp]; %save to one vector
        durClickcon = [durClickcon;durClick];
        nDurcon = [nDurcon; nDur];
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
xlswrite([inDir,'\',choppedDir{3},'_ClicksOnlyConcatCHAR',filedate,'.xls'],clickDnumChar)


%%%Section added to do post-processing where all the clicks are together,
%%%not speparted by xwav. 

%Get detectionTimes
inpath = 'C:\Users\Karlina.Merkens\Documents\Porpoise\OtherRecordings\EJacobson_Harbor_wild\harbor';
infile = 'EJacobsen_harbor_wild_log_150702.xls';
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
specClickTf = specClickTfcon;
specNoiseTf = specNoiseTfcon;
peakFr = peakFrcon;
nDur = nDurcon;
yFilt = yFiltcon;
GraphDir = [inDir,'\matlab_graphs'];


[medianValues,meanSpecClicks,iciEncs] = plotClickEncounters_posthoc_150310(encounterTimes,clickTimes,ppSignal,durClick,...
    specClickTf,specNoiseTf,peakFr,nDur,yFilt,hdr,GraphDir,f);


%Then save everything
save([inDir,'\',choppedDir{3},'_ClicksOnlyConcat',filedate,'.mat'],...
    'clickDnum','durClickcon','nDurcon', 'peakFrcon','ppSignalcon',...
    'specClickTfcon','yFiltcon','medianValues','meanSpecClicks','iciEncs','f')



