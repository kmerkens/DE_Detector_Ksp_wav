function parametersST = dLoad_STsettings

% Assign short term detector settings

parametersST.buff = 500; % # of buffer samples to add on either side of area of interest
parametersST.chan = 1; % which channel do you want to look at?

%parametersST.fRanges = [50000 99000]; 
%parametersST.fRanges = [5000,159000]; % For Kogia on 320
%parametersST.fRanges = [50000,200000]; % For Kogia on 500
%parametersST.fRanges = [50000,187000]; % For Kogia on 375
%parametersST.fRanges = [5000,200000]; % For Dalls on 480
parametersST.fRanges = [5000,143500]; % For Kogia on 288


%parametersST.thresholds = 13500; % Amplitude threshold in counts. 
%parametersST.thresholds = 3500; % Trying to find something that will get 
%the very quiet kogia clicks. 
%150311 - changed from 5000 - was this set wrong? That's higer than the
%3500 in the hi res which it shouldn't be.  are t he unit s the same?
% For predictability, keep this consistent between low and hi res steps.
%parametersST.thresholds = 20000; %Set high for D. Mann wav files, no tf. 
parametersST.thresholds = 50; %Set lower to start V. Janik files
%parametersST.thresholds = 800000; %Set high again for TYack dalls

parametersST.frameLengthSec = .01; %Used for calculating fft size
parametersST.overlap = .50; % fft overlap
parametersST.REWavExt = '(\.x)?\.wav';%  expression to match .wav or .x.wav

% if you're using wav files that have a time stamp in the name, put a
% regular expression for extracting that here:
parametersST.DateRE = '_(\d*)_(\d*)';
% mine look like "filename_20110901_234905.wav" 
% ie "*_yyyymmdd_HHMMSS.wav"


