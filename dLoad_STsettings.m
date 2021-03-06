function parametersST = dLoad_STsettings

% Assign short term detector settings

parametersST.buff = 500; % # of buffer samples to add on either side of area of interest
parametersST.chan = 1; % which channel do you want to look at?

%parametersST.fRanges = [50000 99000]; 
%parametersST.fRanges = [1000,159000]; % For Kogia on 320
% parametersST.fRanges = [10000,170000]; % For Ksima manuscript analysis
% parametersST.fRanges = [10000,227000]; % For Kogia on 500 D Mann
% parametersST.fRanges = [1000,187000]; % For Kogia on 375 V Yanik
% parametersST.fRanges = [1000,190000]; % For Ksima on 384 CARB (CRP)
% parametersST.fRanges = [120000,140000]; % For Ksima on 576 TGridley
%parametersST.fRanges = [5000,200000]; % For Dalls on 480 T Yack
%parametersST.fRanges = [5000,143500]; % For Kogia on 288
%parametersST.fRanges = [5000,249000]; % For Harbor on 500 (E Jacobsen)
% parametersST.fRanges = [90000,143000]; % For DASPR 143
parametersST.fRanges = [1000,15900]; % For C Jenner



%parametersST.thresholds = 13500; % Amplitude threshold in counts. 
%parametersST.thresholds = 3500; % Trying to find something that will get 
%the very quiet kogia clicks. 
%150311 - changed from 5000 - was this set wrong? That's higer than the
%3500 in the hi res which it shouldn't be.  are t he unit s the same?
% For predictability, keep this consistent between low and hi res steps.

% parametersST.thresholds = 40000; %Set high for D. Mann wav files, no tf. 
% parametersST.thresholds = 500; %Set lower to start V. Janik files
% parametersST.thresholds = 2500; %Adjusting for CARB file
% parametersST.thresholds = 20000; %Adjusting for TGridley file
%parametersST.thresholds = 1000000; %Set even higher again for TYack dalls
%parametersST.thresholds = 750000; %Set low to start E Jacobsen porpoise
% parametersST.thresholds = 20; %Adjusting for DASPR
parametersST.thresholds = 2000; %Adjusting for CJenner



parametersST.frameLengthSec = .01; %Used for calculating fft size
parametersST.overlap = .50; % fft overlap
parametersST.REWavExt = '(\.x)?\.wav';%  expression to match .wav or .x.wav
% parametersST.REWavExt = '(\w+).(\.x)?\.wav';%  expression to match .wav or .x.wav

% if you're using wav files that have a time stamp in the name, put a
% regular expression for extracting that here:
parametersST.DateRE = '_(\d*)_(\d*)';
% mine look like "filename_20110901_234905.wav" 
% ie "*_yyyymmdd_HHMMSS.wav"
%Changing to work with DASBR file names, which are
%0123456789.yymmddHHMMSS.wav
% parametersST.DateRE = '.(\d{12})';




