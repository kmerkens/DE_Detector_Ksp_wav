function parametersHR = dLoad_HRsettings

%%% Filter and FFT params %%
%parametersHR.bpRanges = [50000,99000]; % Bandpass filter params in Hz [min,max]
%parametersHR.bpRanges = [5000,159000]; % For Kogia on 320
%parametersHR.bpRanges = [5000,200000]; % For Kogia on 500 (DMann - avoiding noise)
%parametersHR.bpRanges = [1000,187000]; % For Kogia on 375
parametersHR.bpRanges = [5000,200000]; % For dalls on 480 T Yack
%parametersHR.bpRanges = [5000,143500]; %For Kogia on 288

parametersHR.frameLengthUs = 1200; % For fft computation
parametersHR.overlap = .5; % FFT overlap (in decimal, not percent form)
parametersHR.chan = 1; % which channel do you want to look at?
parametersHR.clipThreshold = .98;%  Normalized clipping threshold btwn 0 and 1.  If empty, 
% assumes no clipping.


%%% Recieved level threshold params %%%
%parametersHR.ppThresh = 100;% minimum  RL threshold - dB peak to peak.
parametersHR.ppThresh = 40;% minimum  RL threshold - dB peak to peak. 
%Decreased for DMann wav, particularly because no tf available. 
%parametersHR.ppThresh = 20; %decresed further for V. Janik data.

%parametersHR.countThresh = 3500; % Keep consistent with Lo-res for predictability.
% Can be higher than low res, but not lower!
% Keep count threshold less than equivalent pp threshold. 
%   dBs = 10*log10(abs(fft(counts *2^14))) - 10*log10(fs/(length(fftWindow)))...
%            + transfer function
% note: array uses 2^15
%parametersHR.countThresh = 60000; %Set high for D. Mann wav files. 
%parametersHR.countThresh = 500; %Set lower for V. Janik wav files. 
%parametersHR.countThresh = 1000000; %Set very hi for dalls wav files. 
parametersHR.countThresh = 750000; %set high for EJ Harbor porpoise


%%% Envelope params %%%
parametersHR.energyThr = 0.5; % n-percent energy threshold for envelope duration
parametersHR.dEvLims = [-.4,.9];  % [min,max] Envelope energy distribution comparing 
% first half to second half of high energy envelope of click. If there is
% more energy in the first half of the click (dolphin) dEv >0, If it's more
% in the second half (boats?) dEv<0. If it's about the same (beaked whale)
% dEnv ~= 0 , but still allow a range...
%parametersHR.delphClickDurLims = [5,30];% [min,max] duration in microsec 
% allowed for high energy envelope of click
parametersHR.delphClickDurLims = [5,100];%increased for wav
%parametersHR.delphClickDurLims = [5,100];%lowered min for dalls wav


%%% Other pruning params %%%
%parametersHR.cutPeakBelowKHz = 80; % discard click if peak frequency below X kHz
parametersHR.cutPeakBelowKHz = 100; %% For Kogia on 500 (DMann)
%parametersHR.cutPeakBelowKHz = 110; %% For dalls on 480
%parametersHR.cutPeakAboveKHz = 99.9; % discard click if peak frequency above Y kHz 
parametersHR.cutPeakAboveKHz = 150;%% For Kogia on 500 (DMann)
parametersHR.minClick_us = 16;% Minimum duration of a click in us 
parametersHR.maxClick_us = 1000; % Max duration of a click including echos
%parametersHR.maxNeighbor = 1; % max time in seconds allowed between neighboring 
% clicks. Clicks that are far from neighbors can be rejected using this parameter,
% good for dolphins in noisy environments because lone clicks or pairs of
% clicks are likely false positives
parametersHR.maxNeighbor = 2; %Increased, to get faint kogia detections. 

%parametersHR.mergeThr = 50;% min gap between energy peaks in us. Anything less
% will be merged into one detection the beginning of the next is fewer
% samples than this, the signals will be merged.
parametersHR.mergeThr = 20;%reduced to try to get more shorter clicks in echo bouts. 

% %%%Not necessary for wav data
% parametersHR.localminbottom = 63; %Frequency above which will be checked 
% %for local min
% parametersHR.localmintop = 90; %Frequency below which will be checked for
% %local min.
% parametersHR.localminThr = 68; %Frequency above which a local minimum must
% %exist in order to be thrown out.  If the min is below, the click will be
% %kept


% if you're using wav files that have a time stamp in the name, put a
% regular expression for extracting that here:
parametersHR.DateRE = '_(\d*)_(\d*)';
% mine look like "filename_20110901_234905.wav" 
% ie "*_yyyymmdd_HHMMSS.wav"

%%% Output file extensions. Probably don't need to be changed %%%
parametersHR.clickAnnotExt = 'cTg';
parametersHR.ppExt = 'pTg';
parametersHR.groupAnnotExt = 'gTg';
