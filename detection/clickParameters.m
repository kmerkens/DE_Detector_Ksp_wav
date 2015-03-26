function [clickInd,ppSignal,durClick,bw3db,yNFilt,yFilt,specClickTf,...
    specNoiseTf,peakFr,yFiltBuff,f,deltaEnv,nDur] = clickParameters(~,wideBandData,p,fftWindow,...
    PtfN,clicks,specRange,hdr)

%Take timeseries out of existing file, convert from normalized data to
%counts
%1) calculate spectral received levels RL for click and preceding noise:
%calculate spectra, account for bin width to reach dB re counts^2/Hz,
%add transfer function, compute peak frequency and bandwidth
%2) calculate RLpp at peak frequency: find min & max value of timeseries,
%convert to dB, add transfer function value of peak frequency (should come
%out to be about 9dB lower than value of spectra at peak frequency)
%3) Prune out clicks that don't fall in expected peak frequency, 3dB
%bandwidth/duration range, or which are not high enough amplitude
%(ppSignal)
% ** There's code in here to compute a noise spectrum alongside the click
% spectrum. It should at least get you started if you want that sort of
% thing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables

N = length(fftWindow);
ppSignal = zeros(size(clicks,1),1);
durClick =  zeros(size(clicks,1),1);
bw3db = zeros(size(clicks,1),3);
yNFilt = [];
yFilt = cell(size(clicks,1),1);
zerosvec = zeros(1,300);
yFilt(:) = {zerosvec};
specClickTf = cell(size(clicks,1),1);
yFiltBuff = cell(size(clicks,1),1);
specNoiseTf = cell(size(clicks,1),1);
peakFr = zeros(size(clicks,1),1);
cDLims = ceil([p.minClick_us, p.maxClick_us]./(hdr.fs/1e6));
envDurLim = ceil(p.delphClickDurLims./(hdr.fs/1e6));
nDur = zeros(size(clicks,1),1);
deltaEnv = zeros(size(clicks,1),1);

f = 0:((hdr.fs/2)/1000)/((N/2)-1):((hdr.fs/2)/1000);
f = f(specRange);

% concatonnate vector of noise
% for itr = 1:min([30,size(noiseIn,1)])
%     nStart = noiseIn(itr,1);
%     nEnd = min([noiseIn(itr,2),nStart+580]);
%     yNFilt = [yNFilt,wideBandData(nStart:nEnd)];
% end
% noise = yNFilt*2^14; % convert to counts
buffVal = hdr.fs*.00025; % Add small buffer, here, I want .25 ms, so computing how many samples to use.

validClicks = ones(size(ppSignal));

for c = 1:size(clicks,1)
    % Pull out band passed click timeseries
    yFiltBuff{c} = wideBandData(max(clicks(c,1)-buffVal,1):min(clicks(c,2)+buffVal,size(wideBandData,2)));
    yFiltLength = clicks(c,2)-clicks(c,1)+1;
    yFilt{c,1}(1:yFiltLength) = wideBandData(clicks(c,1):clicks(c,2));
    %convert  timeseries into counts
    if strcmp(hdr.fType, 'xwav')
        click = yFilt{c}*2^14;
        clickBuff = yFiltBuff{c}*2^14;
    else
        click = yFilt{c}*2^15; % array needs 2^15 to convert into counts. 
        % Unclear what the cause of this is.
        clickBuff = yFiltBuff{c}*2^15;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate duration in seconds
    durClick(c) = (clicks(c,2)-clicks(c,1));
    
    % Compute click spectrum
    winLength = length(clickBuff);
    wind = hann(winLength);
    wClick = zeros(1,N);
    wClick(1:winLength) = clickBuff.*wind.';
    spClick = 10*log10(abs(fft(wClick,N)));
    
    % Compute noise spectrum
%     windNoise = hann(N);
%     wNoise = [];
%     for itr1 = 1:1:floor(length(noise)/N)
%         wNoise(:,itr1) = noise(1,((itr1-1)*N)+1:((itr1-1)*N)+N).*windNoise.';
%     end
%     spNoise = fastsmooth(mean(10*log10(abs(fft(wNoise,N))),2),15,1,1)';
    
    % account for bin width
    sub = 10*log10(hdr.fs/N);
    spClickSub = spClick-sub;
    % spNoiseSub = spNoise-sub;
    
    %reduce data to first half of spectra
    spClickSub = spClickSub(:,1:N/2);
    % spNoiseSub = spNoiseSub(:,1:N/2);
    
    specClickTf{c} = spClickSub(specRange)'+PtfN;
    % specNoiseTf{c} = spNoiseSub(specRange)'+PtfN;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % calculate peak click frequency
    % max value in the first half samples of the spectrogram
    
    [valMx, posMx] = max(specClickTf{c}(1:end-1)); % ignore last point - bandpass issue
    peakFr(c) = f(posMx); %peak frequency in kHz
    
    if peakFr(c) < p.cutPeakBelowKHz || peakFr(c) > p.cutPeakAboveKHz
        validClicks(c) = 0;
        continue
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate RLpp at peak frequency: find min/max value of timeseries,
    %convert to dB, add transfer function value of peak frequency (should come
    %out to be about 9dB lower than value of spectra at peak frequency)
    
    % find lowest and highest number in timeseries (counts) and add those
    high = max(yFilt{c}.');
    low = min(yFilt{c}.');
    ppCount = high+abs(low);
    
    %calculate dB value of counts and add transfer function value at peak
    %frequency to get ppSignal (dB re 1uPa)
    P = 20*log10(ppCount);
    
    peakLow=floor(peakFr(c));
    if peakLow == (hdr.fs/2)/1000
        fLow = N/2;
    else
        fLow=find(f>peakLow);
    end
    
    %add PtfN transfer function at peak frequency to P
    tfPeak = PtfN(fLow(1));
    ppSignal(c) = P+tfPeak;
    
    if ppSignal(c)< p.ppThresh
        validClicks(c) = 0;
        continue
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Not necessary for 320Khz Data
%     %Use the ratio of the median energy near 75 and 97kHz as a cutoff,
%     %throwing out clicks with too low a ratio (that means there's too much
%     %energy at lower frequencies/likely a click lower down)
%     lowsetC = find(f>72,1,'first');
%     lowsetD = find(f<77,1,'last');
%     highsetC = find(f>94,1,'first');
%     highsetD = find(f<99,1,'last');
%     %Get the median from each
%     lowmed = median(specClickTf{c}(lowsetC:lowsetD));
%     highmed = median(specClickTf{c}(highsetC:highsetD));
%     specratioH(c) = highmed/lowmed; 
% 
%     if specratioH(c) < 1.07
%         validClicks(c) = 0;
%         continue
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Not necessary for 320Khz Data
%     %Use the ratio of the median energy near 55 kHz and 70 kHz as a cutoff
%     %- if it's lower than 1.15, throw the click out. 
%     %identify the bins that will be used to find the low and high values
%     %for ratio calculation
%     lowsetA = find(f>52,1,'first');
%     lowsetB = find(f<57,1,'last');
%     highsetA = find(f>67,1,'first');
%     highsetB = find(f<72,1,'last');
%     %Get the median from each
%     lowmed = median(specClickTf{c}(lowsetA:lowsetB));
%     highmed = median(specClickTf{c}(highsetA:highsetB));
%     specratioL(c) = highmed/lowmed; 
%     
%     if specratioL(c) < 1.165
%         validClicks(c) = 0;
%         continue
%     end
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%Not necessary for 320Khz Data
%     %check for local minimum between 60 and 90 kHz, and if it exists closer
%     %to 90, the click is not valid (likely dsp).
%     %Find the bins that match the frequencies 60-80 kHz
%     bottombin = find(f >= p.localminbottom,1, 'first');
%     topbin = find(f <p.localmintop,1, 'last');
%     %Find the minimum in that range
%     localmindb = min(specClickTf{c}(bottombin:topbin));
%     localminbin = find((specClickTf{c}(1:end-1)) == localmindb);
%     localminf(c) = f(localminbin);
% 
%     if localminf(c) > p.localminThr
%         validClicks(c) = 0;
%         continue
%     end
    
   
    %%%%%%%%%%%%%%%%%
    % calculate click envelope (code & concept from SBP 2014):
    % pre_env = hilbert(click);
    % env = sqrt((real(pre_env)).^2+(imag(pre_env)).^2); %Au 1993, S.178, equation 9-4
    env = abs(hilbert(click));
    
    %calculate energy duration over x% energy
    env = env - min(env);
    env = env/max(env);
    %determine if the slope of the envelope is positive or negative
    %above x% energy
    aboveThr = find(env>=p.energyThr);
    direction = nan(1,length(aboveThr));
    for a = 1:length(aboveThr)
        if aboveThr(a)>1 && aboveThr(a)<length(env)
            % if it's not the first or last element fo the envelope, then
            % -1 is for negative slope, +1 is for + slope
            delta = env(aboveThr(a)+1)-env(aboveThr(a));
            if delta>=0
                direction(a) = 1;
            else
                direction(a) = -1;
            end
        elseif aboveThr(a) == 1
            % if you're looking at the first element of the envelope
            % above the energy threshold, consider slope to be negative
            direction(a) = -1;
        else  % if you're looking at the last element of the envelope
            % above the energy threshold, consider slope to be positive
            direction(a) = 1;
        end
    end
    
    %find the first value above threshold with positive slope and find
    %the last above with negative slope
    lowIdx = aboveThr(find(direction,1,'first'));
    negative = find(direction==-1);
    if isempty(negative)
        highIdx = aboveThr(end);
    else
        highIdx = aboveThr(negative(end));
    end
    nDur(c,1) = highIdx - lowIdx + 1;
    
    if nDur(c)>(p.delphClickDurLims(2))
        validClicks(c) = 0;
        continue
    end
    
    %compare maximum first half of points with second half.
    halves = ceil(nDur(c,1)/2);   
    env1max = max(env(lowIdx:min([lowIdx+halves,length(env)])));
    env2max = max(env(min([lowIdx+(halves)+1,length(env)]):end));
    deltaEnv(c,1) = env1max-env2max;
    
     if deltaEnv(c) < p.dEvLims(1)
        validClicks(c) = 0;
        continue
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate duration
    durClick(c) = clicks(c,2)-clicks(c,1);
    
    if durClick(c) > cDLims(2)
        validClicks(c) = 0;
        continue
    end
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Make a plot of the spectrum, to show what we're dealing with
%     figure
%     plot(f(1:end-1),specClickTf{c}(1:end-1))
%     close

    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %calculate bandwidth
%     %-3dB bandwidth
%     %calculation of -3dB bandwidth - amplitude associated with the halfpower points of a pressure pulse (see Au 1993, p.118);
%     low = valMx-3; %p1/2power = 10log(p^2max/2) = 20log(pmax)-3dB = 0.707*pmax; 1/10^(3/20)=0.707
%     %walk along spectrogram until low is reached on either side
%     slopeup=fliplr(specClickTf{c}(1:posMx));
%     slopedown=specClickTf{c}(posMx:round(length(specClickTf{c})));
%     for e3dB=1:length(slopeup)
%         if slopeup(e3dB)<low %stop at value < -3dB: point of lowest frequency
%             break
%         end
%     end
%     for o3dB=1:length(slopedown)
%         if slopedown(o3dB)<low %stop at value < -3dB: point of highest frequency
%             break
%         end
%     end
%     
%     %calculation from spectrogram -> from 0 to 100kHz in 256 steps (FFT=512)
%     high3dB = (hdr.fs/(2*1000))*((posMx+o3dB)/(length(specClickTf{c}))); %-3dB highest frequency in kHz
%     low3dB = (hdr.fs/(2*1000))*(posMx-e3dB)/(length(specClickTf{c})); %-3dB lowest frequency in kHz
%     bw3 = high3dB-low3dB;
%     
%     bw3db(c,:)= [low3dB, high3dB, bw3];
%     
%     
%     %%%%%
    
    
end



% validClicks = ones(size(ppSignal));
% 
% % Check parameter values for each click
% for idx = 1:length(ppSignal)
%     tfVec = [deltaEnv(idx) < p.dEvLims(1);...
%              peakFr(idx) < p.cutPeakBelowKHz;...
%              peakFr(idx) > p.cutPeakAboveKHz;...
%              nDur(idx)>  (envDurLim(2));...
%              nDur(idx)<  (envDurLim(1));...
%              durClick(idx) > cDLims(2)];
% %          plot(yFiltBuff{idx})
% %          title(sum(tfVec))
% %          1;
%     if ppSignal(idx)< p.ppThresh
%         validClicks(idx) = 0; 
%     elseif sum(tfVec)>0   
%         validClicks(idx) = 0; 
%     end
%     
% end
clickInd = find(validClicks == 1);
% throw out clicks that don't fall in desired ranges
ppSignal = ppSignal(clickInd,:);
durClick =  durClick(clickInd,:);
% bw3db = bw3db(clickInd,:);
% frames = frames{clickInd};
% yFilt = {yFilt{clickInd}};
yFilt = yFilt(clickInd);
yFiltBuff = yFiltBuff(clickInd);
specClickTf = specClickTf(clickInd);
% specNoiseTf = {specNoiseTf{clickInd}};
peakFr = peakFr(clickInd,:);
% yNFilt = {yNFilt};
deltaEnv = deltaEnv(clickInd,:);
nDur = nDur(clickInd,:);
