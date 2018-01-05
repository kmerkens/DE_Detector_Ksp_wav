function [clickInd,ppSignal,durClick,dur95usec,dur95usecTails,durRMSus,bw3db,...
    bw10db,bwRMS,QRMS,Q3dB,yNFilt,yFilt,specClickTf,...
    specNoiseTf,peakFr,centFr,snr,yFiltBuff,f,deltaEnv,...
    nDur] = clickParameters(noises,wideBandData,p,fftWindow,...
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
dur95usec = zeros(size(clicks,1),1);
dur95usecTails = zeros(size(clicks,1),1);
durRMSus = zeros(size(clicks,1),1);
bw3db = zeros(size(clicks,1),3);
bw10db = zeros(size(clicks,1),3);
bwRMS = zeros(size(clicks,1),1);
QRMS = zeros(size(clicks,1),1);
Q3dB = zeros(size(clicks,1),1);
yNFilt = cell(size(clicks,1),1);
yFilt = cell(size(clicks,1),1);
% zerosvec = zeros(1,300);
% yFilt(:) = {zerosvec};
specClickTf = cell(size(clicks,1),1);
yFiltBuff = cell(size(clicks,1),1);
specNoiseTf = cell(size(clicks,1),1);
yNFiltBuff = cell(size(clicks,1),1);
peakFr = zeros(size(clicks,1),1);
centFr = zeros(size(clicks,1),1);
snr = zeros(size(clicks,1),1);
cDLims = ceil([p.minClick_us, p.maxClick_us]./(hdr.fs/1e6));
cDLims95 = ceil([p.minClick_us, p.maxClick95_us]./(hdr.fs/1e6));
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

% buffVal = hdr.fs*.00025; % Add small buffer, here, I want 250 us ms, so computing how many samples to use.
%For calculating duration
buffValBefore = hdr.fs*.00001; %adding 10us buffer before
buffValAfter = hdr.fs*.00002; %adding 20us buffer after
%buffVal = buffVal * 4;
validClicks = ones(size(ppSignal));


for c = 1:size(clicks,1)
    % Pull out band passed click timeseries
%     yFiltBuff{c} = wideBandData(floor(max(clicks(c,1)-buffVal,1)):ceil(min(clicks(c,2)+buffVal,size(wideBandData,2))));
%     yNFiltBuff{c} = wideBandData(floor(max(noises(c,1))):ceil(min(noises(c,2)+buffVal,size(wideBandData,2))));
    yFiltBuff{c} = wideBandData(floor(max(clicks(c,1)-buffValBefore,1)):ceil(min(clicks(c,2)+buffValAfter,size(wideBandData,2))));
    yNFiltBuff{c} = wideBandData(floor(max(noises(c,1))):ceil(min(noises(c,2)+buffValAfter,size(wideBandData,2))));
    yFiltLength = clicks(c,2)-clicks(c,1)+1;
    yNFiltLength = noises(c,2)-noises(c,1)+1;
%     yFilt{c,1}(1:yFiltLength) = wideBandData(clicks(c,1):clicks(c,2)); %changed to only extract portion of bandwidth allowed through BPF
%     yNFilt{c,1}(1:yNFiltLength) = wideBandData(noises(c,1):noises(c,2));
    yFilt{c} =  wideBandData(clicks(c,1):clicks(c,2));
    yNFilt{c} = wideBandData(noises(c,1):noises(c,2));
    %convert  timeseries into counts
    if strcmp(hdr.fType, 'xwav')
        click = yFilt{c}*2^14;
        clickBuff = yFiltBuff{c}*2^14;
        noise = yNFilt{c}*2^14;
        noiseBuff = yNFiltBuff{c}*2^14;
    else
        click = yFilt{c}*2^15; % array needs 2^15 to convert into counts. 
        % Unclear what the cause of this is.
        clickBuff = yFiltBuff{c}*2^15;
        noise = yNFilt{c}*2^15;
        noiseBuff = yNFiltBuff{c}*2^15;
    end
        
    % Compute click spectrum
    winLength = length(clickBuff);
    wind = hann(winLength);
    %wClick = zeros(1,N);
    wClick(1:winLength) = clickBuff.*wind';
    %wClick(1:winLength) = clickBuff.*wind.'; %KPM for two vectors it seems
    %that whether the second . is there or not does not matter, because the 
    %type of transpose doesn't matter. 
    %wClick(1:winLength) = click.*wind.';
    spClick = 20*log10(abs(fft(wClick,N)));
    
    % Compute noise spectrum
    %%%karlis
    if length(clickBuff) < length(noiseBuff);
        noiseBuff = noiseBuff(1:length(clickBuff));%make it the same length as the click
    else
        %display('Uh oh! The lengths of the click and noise are not right!')
        validClicks(c) = 0; %remove this click/noise pair from analysis
        continue
    end
    %winNLength = length(noiseBuff);
    winNLength = length(clickBuff); %make it the same length as the click
    %winNLength = length(noise);
    windN = hann(winNLength);
    %wNoise = zeros(1,N);
    wNoise(1:winNLength) = noiseBuff.*windN';
    %wNoise(1:winNLength) = noiseBuff.*windN.';%KPM I don't think that second .
    %should be there after wind
    %wNoise(1:winNLength) = noise.*windN.';
    spNoise = 20*log10(abs(fft(wNoise,N)));
    %spNoise = fastsmooth(20*log10(abs(fft(wNoise,N))),15,1,1); %smooth is prettier
    
%     %%%kaits
%     windNoise = hann(N);
%     wNoise = [];
%     for itr1 = 1:1:floor(length(noise)/N)
%         wNoise(:,itr1) = noise(1,((itr1-1)*N)+1:((itr1-1)*N)+N).*windNoise.';
%     end
%     spNoise = fastsmooth(mean(20*log10(abs(fft(wNoise,N))),2),15,1,1)';
%     
    % account for bin width
    sub = 10*log10(hdr.fs/N);
    spClickSub = spClick-sub;
    spNoiseSub = spNoise-sub;
    
    %reduce data to first half of spectra
    spClickSub = spClickSub(:,1:N/2);
    spNoiseSub = spNoiseSub(:,1:N/2);
    
    %Add TF
    specClickTf{c} = spClickSub(specRange)'+PtfN;
    specNoiseTf{c} = spNoiseSub(specRange)'+PtfN;
    
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Look at plot of click to see where it falls in the timeseries and then
%     %follow why it fails. 
%      subplot(2,1,1)
%      plot(wideBandData)
%      hold on
%      line([clicks(c,1) clicks(c,1)], [min(wideBandData) max(wideBandData)],'Color','r','LineWidth',1);
%      line([clicks(c,2) clicks(c,2)], [min(wideBandData) max(wideBandData)],'Color','r','LineWidth',1);
%      
%      subplot(2,1,2)
%      startshortplot = clicks(c,1)-100;
%      endshortplot = clicks(c,2)+100;
%      plotseg = wideBandData(1,startshortplot:endshortplot);
%      plot(plotseg);
%      line([100 100], [min(wideBandData) max(wideBandData)],'Color','r','LineWidth',1);
%      line([size(plotseg,2)-100 size(plotseg,2)-100], [min(wideBandData) max(wideBandData)],'Color','r','LineWidth',1);
%      
%      close(figure(1))
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % calculate peak click frequency
    % max value in the first half samples of the spectrogram
    
    [valMx, posMx] = max(specClickTf{c}(1:end-1)); % ignore last point - bandpass issue
    peakFreq = f(posMx); %peak frequency in kHz
    peakFr(c) = peakFreq;
    
    if peakFr(c) < p.cutPeakBelowKHz || peakFr(c) > p.cutPeakAboveKHz
        validClicks(c) = 0;
        continue
    end
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %if peakFr(c) > 140; %make a plot of the spectrum to check what's up 
%     %with these really high freq signals
%          plot(f(1:end-1),specClickTf{c}(1:end-1))
%     %end
        
    
    
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate the centroid frequency of the click, using Simone's code
    %frequency centroid (or center frequency) in kHz. Au 1993 Equation 10-3
    linearSpec = 10.^(specClickTf{c}/20)'; %undo the 20*log that was done earlier
    %Freq_vec is the frequency band of the spectrum, with steps according
    %to sample rate, which is the same thing as f, calculated above.
%     Freq_vec=0:(hdr.fs/2)/(length(linearSpec)-1):0.5*hdr.fs; %The vector of frequencies = m
    Freq_vec = f*1000; %The correct frequency bins, in Hz.
    centroidFr = (sum(Freq_vec.*(linearSpec.^2))/sum(linearSpec.^2))/1000; %Calculate the cf, /1000 to get to kHz
    centFr(c) = centroidFr;
    %For now, keep all centFr values - don't use this for excluding any
    %clicks. 
    
%     if peakFr < 120
%     	figure(1)
%         plot(f(1:end-1),linearSpec);
%         figure(2)
%         plot(f(1:end-1),specClickTf{c}(1:end-1));
%     end
%     
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate the RMS bandwidth in kHz, using Simone's code
    %RMS bandwidth in kHz
    nc_msbandw=(sum(Freq_vec.^2.*(linearSpec.^2))/sum(linearSpec.^2)); % Au 1993, equation 10-4, 
    centroidFrHZ = centroidFr*1000; %in Hz
    ms_bandw=nc_msbandw-centroidFrHZ.^2; % Au 1993, equation 10-6
    %nc_rmsbandw=sqrt(nc_msbandw);
    rmsBW=(sqrt(ms_bandw))/1000; % Root of mean square bandwidth -> RMS /1000 to get to kHz
    
    bwRMS(c) = rmsBW;
%     figure(1)
%     plot(linearSpec)
%     figure(2)
%     plot(specClickTf{c})
%     
    

    
    
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %From Simone's code:
%     %calculate further yFilt parameters
    
%     %frequency centroid (or center frequency) in kHz
%     linearSpec=10.^(click/20); %click = the first half of the spectra = specClickTf{c}, this undoes the 20 * log 10 that was implemented earlier
%     Freq_vec=0:(fs/2)/(length(linearSpec)-1):0.5*fs; %The vector of frequencies = m
%     F0(n)=(sum(Freq_vec.*linearSpec(1:length(linearSpec)).^2)/sum(linearSpec(1:length(linearSpec)).^2))/1000; % Au 1993, equation 10-3
%     %KPM I'm not sure what the ^2 and the /1000 mean, but I trust the code, and
%     %the equation they cite, and the output is reasonable.

%     %RMS bandwidth in kHz
%     nc_msbandw=(sum(Freq_vec.^2.*linearSpec(1:length(linearSpec)).^2)/sum(linearSpec(1:length(linearSpec)).^2)); % Au 1993, equation 10-4
%     ms_bandw=nc_msbandw-F0(n).^2; % Au 1993, equation 10-6 %n is the click number for the loop, so my {c}, F0 is the centroid frequency, calculated above
%     %nc_rmsbandw=sqrt(nc_msbandw);
%     rmsBandw(n)=(sqrt(ms_bandw))/1000; % Root of mean square bandwidth -> RMS

%     %time centroid in s = centroid of the time waveform
%     yFiltDur=dur(n)/1000;
%     tClick=floor(yFiltDur*fs);
%     t=0:yFiltDur/(tClick-1):yFiltDur;
%     t0(n)=sum(t.*tClick.^2)/sum(tClick.^2); % Au 1993, equation 10-29
% 
%     %rms duration in s
%     rmsDuration(n)=sqrt(sum((t-t0(n)).^2.*tClick.^2)/sum(tClick.^2));



    %%%%%%%%%%%%%%%%%
    % calculate click envelope (code & concept from SBP 2014):
    % pre_env = hilbert(click);
    % env = sqrt((real(pre_env)).^2+(imag(pre_env)).^2); %Au 1993, S.178, equation 9-4
    env1 = abs(hilbert(click));
    %env1 = abs(hilbert(clickBuff)); %Not using with buff, because this is
    %also used for calculating duration, and we don't want that buffer in
    %the duration calculation.
    
    %calculate energy duration over x% energy (normalize to 0:1)
    env2 = env1 - min(env1);
    env3 = env2/max(env2);
    %determine if the slope of the envelope is positive or negative
    %above x% energy
    aboveThr = find(env3>=p.energyThr);
    direction = nan(1,length(aboveThr));
    for a = 1:length(aboveThr)
        if aboveThr(a)>1 && aboveThr(a)<length(env3)
            % if it's not the first or last element fo the envelope, then
            % -1 is for negative slope, +1 is for + slope
            delta = env3(aboveThr(a)+1)-env3(aboveThr(a));
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
    
    %%%%%%%%%%%%%%%%%
    % %find the first value above threshold with positive slope and find
    % %the last above with negative slope
    % lowIdx = aboveThr(find(direction,1,'first'));
    % negative = find(direction==-1);
    % if isempty(negative)
        % highIdx = aboveThr(end);
    % else
        % highIdx = aboveThr(negative(end));
    % end
    % nDur(c,1) = highIdx - lowIdx + 1;
    
    % if nDur(c)>(p.delphClickDurLims(2))
        % validClicks(c) = 0;
        % continue
    % end
    
    % %compare maximum first half of points with second half.
    % halves = ceil(nDur(c,1)/2);   
    % env1stmax = max(env3(lowIdx:min([lowIdx+halves,length(env3)])));
    % env2ndmax = max(env3(min([lowIdx+(halves)+1,length(env3)]):end));
    % deltaEnv(c,1) = env1stmax-env2ndmax;
    
     % if deltaEnv(c) < p.dEvLims(1)
        % validClicks(c) = 0;
        % continue
     % end
     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %New section to calculate the duration of the click using 97% of of
     %the envelope as the metric, sensu madsen et al 2005 and 2004,
     %etc. That means looking from 1.5% to 98.5%

     %Try making cumulative energy plot. 
     %Standardize out of the total energy
     
     %Then sum and normalize to the max  
     totenvenergy = sum(env1);
     env4 = env1/totenvenergy;
     clear('envcum')
     for e = 2:size(env4,2)
        envcum(1) = 0;
        envcum(e) = envcum(e-1)+env4(e-1);
        if e == size(env4,2)
            envcum(e+1) = envcum(e) + env4(e);
        end
     end
%      %Locate the ends of the 97% range
%      durStart = 0.015;
%      durEnd = 0.985;
%      afterStart = find(envcum<=durStart,1,'first'); 
%      beforeEnd = find(envcum>=durEnd,1,'first');
%      dur97cts = beforeEnd - afterStart+1;
%      dur97usec1 = dur97cts/(hdr.fs/1e6); %just to be able to see it while testing
%      dur97usec(c) = dur97cts/(hdr.fs/1e6);
     
     %%Try for 95%
     durStart95 = 0.025;
     durEnd95 = 0.975;
     afterStart95 = find(envcum<=durStart95,1,'first'); 
     beforeEnd95 = find(envcum>=durEnd95,1,'first');
     dur95cts = beforeEnd95 - afterStart95+1;
     
     dur95usec1 = dur95cts/(hdr.fs/1e6); %just to be able to see it while testing
     dur95usec(c) = dur95cts/(hdr.fs/1e6);
     %Only save durations less than 500 us, or whatever is specified in the
     %HR settings
     if dur95usec(c) > p.maxClick95_us
         validClicks(c) = 0;
         continue
     end
        
     %Test going from each end to find the 1.5 and 98.5% points. I'll call
     %this the "tails" version
     durTails95 = 0.025;
     afterStartTails95 = find(env3>=durTails95,1,'first'); 
     beforeEndTails95 = find(env3>=durTails95,1,'last');
     dur95ctsTails = beforeEndTails95 - afterStartTails95+1;
     dur95usecTails1 = dur95ctsTails/(hdr.fs/1e6);
     dur95usecTails(c) = dur95ctsTails/(hdr.fs/1e6);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Calculate duration
    
     durClick(c) = clicks(c,2)-clicks(c,1);
     clickdur = durClick(c)/(hdr.fs/1e6);
     if durClick(c) > cDLims(2)
        validClicks(c) = 0;
        continue
     end

%      %Check to see how much the extra buffer might be influencing the total
%      %duration of the click
%      figure(1)
%      subplot(1,2,1)
% %      plot(clickBuff)
%      plot(click)
%      hold on
%      plot(env1,'r')
% %      line([afterStart95 afterStart95], [min(clickBuff) max(clickBuff)],'Color','k','LineWidth',1);
% %      line([beforeEnd95 beforeEnd95], [min(clickBuff) max(clickBuff)],'Color','k','LineWidth',1);
% %      text(1,min(clickBuff),['dur95cts=',num2str(dur95cts)]);
% %      text(1,min(clickBuff)-8e5,['dur95us=',num2str(dur95usec1)]);
%      line([afterStart95 afterStart95], [min(click) max(click)],'Color','k','LineWidth',1);
%      line([beforeEnd95 beforeEnd95], [min(click) max(click)],'Color','k','LineWidth',1);
%      text(1,min(click),['dur95cts=',num2str(dur95cts)]);
%      text(1,min(click)-8e5,['dur95us=',num2str(dur95usec1)]);     
%      
% %      %Comparing with durClick
% %      subplot(2,1,2)
% %      startshortplot = clicks(c,1)-1;
% %      endshortplot = clicks(c,2)+10;
% %      plotseg = wideBandData(1,startshortplot:endshortplot);
% %      plot(plotseg);
% %      line([1, 1], [min(plotseg) max(plotseg)],'Color','r','LineWidth',1);
% %      line([size(plotseg,2)-10 size(plotseg,2)-10], [min(plotseg) max(plotseg)],'Color','r','LineWidth',1);
% %      %text(1,max(wideBandData)+20,['durClick-cts=',num2str(durClick(c))]);
% %      text(1,max(wideBandData),['durClick-us=',num2str(clickdur)]);
% %      text(size(plotseg,2)/2,2*max(wideBandData),['difference btwn=',num2str(dur95usec1-clickdur),' us']);
%      
%      %Look at wider window
%      subplot(1,2,2)
%      startshortplot = clicks(c,1)-200;
%      endshortplot = clicks(c,1)+400;
%      plotseg = wideBandData(1,startshortplot:endshortplot);
%      plot(plotseg(1:600));
%      line([200, 200], [min(plotseg) max(plotseg)],'Color','r','LineWidth',1);
%      line([durClick(c)+200, durClick(c)+200], [min(plotseg) max(plotseg)],'Color','r','LineWidth',1);
     
%      %Other options:
%      subplot(2,1,2)
%      plot(click(1,1:150))
%      
%      %check to see where click is within the segment
%      figure(2)
%      subplot(2,1,1)
%      plot(wideBandData)
%      hold on
%      line([clicks(c,1) clicks(c,1)], [min(wideBandData) max(wideBandData)],'Color','r','LineWidth',1);
%      line([clicks(c,2) clicks(c,2)], [min(wideBandData) max(wideBandData)],'Color','r','LineWidth',1);
%      
%      close(figure(1))
%      close(figure(2))
%        


%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Calculate the RMS duration in s, using Simone's code
%     %time centroid in s = centroid of the time waveform
%     %This doesn't really work since the window size of each click varies.
%     %But we can apply the methods below.
%     durClickBuff = size(yFiltBuff{c},2);
%     clickDur = durClickBuff/(hdr.fs/1000); %trying with clickBuff, slightly longer
% %     clickDur = durClick(c)/(hdr.fs/1000); %total click duration in seconds. 
%     tClick = click; %The time series of the waveform, click, in counts
%     clickLength_samp = size(tClick,2);
%     
%     t = 0:clickDur/(clickLength_samp-1):clickDur; %Dividing the click into one time bin per sample
%     timecent = sum(t.*tClick.^2)/sum(tClick.^2); % Au 1993, equation 10-29, time centroid in seconds
% 
%     %rms duration in s
%     rmsDuration(c) = sqrt(sum((t-timecent).^2.*tClick.^2)/sum(tClick.^2)); 
%     rmsDur_ms = rmsDuration(c)*1000;
%     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate the duration in s using a fixed 500 us window around the
    %peak of the timeseries. Hybrid methods from Madsen et al 2005a and 
    %madsen and wahlberg 2007, using sbp code.
    
    %Start by finding peak of timeseries - use the envelope around click to
    %find the peak value and it's location in click
    [tsPeakVal, tsPeakLoc] = max(env1);
    %Use the location, in samples relative to the start of click, to
    %identify the location, in samples relative to the start of this whole
    %WideBandData segment. clicks(c,1) is the start of this click
    tsPeakTime_WBD = clicks(c,1)+ tsPeakLoc;
    %Calculate the number of samples in 500us
    durWinLen = floor(300*hdr.fs/1e6); %192 samples, half is 96
    %Get a new click time series, going 96 samples in each direction from
    %the peak
    durWinStart = ceil(tsPeakTime_WBD - (durWinLen/2));
    durWinEnd = floor(tsPeakTime_WBD + (durWinLen/2));
    durWinBounds = [durWinStart, durWinEnd];
    if durWinEnd > size(wideBandData,2)
        durWinEnd = size(wideBandData,2);
    end
    durWinTS = wideBandData(durWinStart:durWinEnd-1);
    %Use this time series and generate the envelope. 
    env5 = abs(hilbert(durWinTS)); %This is the analytic signal
    %Integrate to get the area under the curve. cumtrapz returns a vector
    %with the area under the 192 segments
    env5Int = cumtrapz(env5);
    env5max = max(env5Int); %the total energy
    %Calculate the values of 2.5% and 97.5% of the total, as the bounds on
    %the 97% total energy
    low025perc = 0.025 * env5max;
    high975perc = 0.975 * env5max;
    %Find the bins wherein this % lies, and round toward the centeroid.
    startbin = find(low025perc > env5Int,1,'last'); 
    endbin = find(high975perc > env5Int,1,'last');
    %Calculate the time between these two bins
    dur500Win_us = (endbin-startbin)/hdr.fs*1e6;
    
    %Test plots to check calculations are below - use the snr to only look
    %at "good clicks"
    
    
    
    %%%%%Calculate the rms duration according to Au (using the 500 us
    %%%%%window from above)
    %Now use this time series to calculate the rms duration. 
    %durWinTS is equivalent to yFilt, so first convert to counts to make something
    %equivalent to click. 
    clickTSdur = durWinTS*2^15;
    %Convert the window duration into seconds
    clickDur = durWinLen/(hdr.fs); %trying with clickBuff, slightly longer
%     clickDur = durClick(c)/(hdr.fs/1000); %total click duration in seconds. 
    tClick = clickTSdur; %The time series of the waveform, click, in counts
    clickLength_samp = size(tClick,2); %This shuold be 192, because we set it above
    t = 0:clickDur/(clickLength_samp-1):clickDur; %Dividing the click into one time bin per sample
    timecent = sum(t.*tClick.^2)/sum(tClick.^2); % Au 1993, equation 10-29, time centroid in seconds
    %should be right around 250 because we placed the window over the
    %signal peak

    %rms duration in s
    durRMSs = sqrt(sum((t-timecent).^2.*tClick.^2)/sum(tClick.^2)); 
    durRMSus(c) = durRMSs*1e6;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Something to assess "shape" of time series - check for "teardrop"
%     %kogia shape - if the number of samples after the peak is larger than
%     %from before. Two thresholds to set before: the bigenv (big envelope),
%     %and the propthresh 
%     %ONLY USE when you want to keep the BEST, PRETTIEST, CLOSEST TO ON AXIS
%     %clicks, like when you want an "ideal" data set, or are comparing
%     %with other, similar species, like Dall's
%     
%     bigenv = 0.01;
%     overthresh = find(env>0.1);
%     numover = size(overthresh,2);
%     smaller = overthresh<(find(env==1));
%     numsmaller = sum(smaller);
%     propsmaller(c) = numsmaller/numover;
%     propthresh = 0.4;
%     
%      if propsmaller(c) > propthresh
%         validClicks(c) = 0;
%         continue
%     end
%     
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate bandwidth
    %-3dB bandwidth
    %calculation of -3dB bandwidth - amplitude associated with the halfpower points of a pressure pulse (see Au 1993, p.118);
    low = valMx-3; %p1/2power = 10log(p^2max/2) = 20log(pmax)-3dB = 0.707*pmax; 1/10^(3/20)=0.707
    %walk along spectrogram until low is reached on either side
%    slopeup=fliplr(specClickTf(c,1:posMx)); Use this if the spectra are
%    %matricies, otherwise use the one below, with flipud, if they're
%    %vectors in cell arrays. 
    slopeup=flipud(specClickTf{c}(1:posMx));

    slopedown=specClickTf{c}(posMx:floor(length(specClickTf{c})));
    for e3dB=1:length(slopeup)
        if slopeup(e3dB)<low %stop at value < -3dB: point of lowest frequency
            break
        end
    end
    for o3dB=1:length(slopedown)
        if slopedown(o3dB)<low %stop at value < -3dB: point of highest frequency
            break
        end
    end
    
%     %calculation from spectrogram -> from 0 to 100kHz in 256 steps (FFT=512)
%     high3dB = (hdr.fs/(2*1000))*((posMx+o3dB)/(length(specClickTf{c}))); %-3dB highest frequency in kHz
%     low3dB = (hdr.fs/(2*1000))*(posMx-e3dB)/(length(specClickTf{c})); %-3dB lowest frequency in kHz
%     bw3 = high3dB-low3dB;
    
    %Alternative method - use f vector previously determined:
    highbin = posMx+o3dB;
    %Add something to remove values that are above the nyquist. bw3db is a
    %vector of NaNs to start with, so in this case we just skip the part of
    %saving those numbers and leave the NaNs.
    if highbin > length(f)
        continue
    else
        high3dB = f(highbin);
        lowbin = posMx-e3dB;
        if lowbin==0
            lowbin=1; %change to lowest vaue possible
        end
        low3dB = f(lowbin);
        bw3 = high3dB-low3dB;
        
        bw3db(c,:)= [low3dB, high3dB, bw3]; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     %calculate bandwidth
    %-10dB bandwidth
    %calculation of -10dB bandwidth - amplitude associated with the halfpower points of a pressure pulse (see Au 1993, p.118);
    low10 = valMx-10; %p1/2power = 10log(p^2max/2) = 20log(pmax)-10dB = 0.3162*pmax; 1/10^(10/20)=0.3162
    %walk along spectrogram until low is reached on either side
    %    slopeup=fliplr(specClickTf(c,1:posMx)); Use this if the spectra are
%    %matricies, otherwise use the one below, with flipud, if they're
%    %vectors (in cell arrays). 
    slopeup10=flipud(specClickTf{c}(1:posMx));

    slopedown10=specClickTf{c}(posMx:floor(length(specClickTf{c})));
    for e10dB=1:length(slopeup10)
        if slopeup10(e10dB)<low10 %stop at value < -10dB: point of lowest frequency
            break
        end
    end
    for o10dB=1:length(slopedown10)
        if slopedown10(o10dB)<low10 %stop at value < -10dB: point of highest frequency
            break
        end
    end
    
%     %calculation from spectrogram -> from 0 to 100kHz in 256 steps (FFT=512)
%     high10dB = (hdr.fs/(2*1000))*((posMx+e10dB)/(length(specClickTf{c}))); %-3dB highest frequency in kHz
%     low10dB = (hdr.fs/(2*1000))*(posMx-e10dB)/(length(specClickTf{c})); %-3dB lowest frequency in kHz
%     bw10 = high10dB-low10dB;
    
    %Alternative method - use f vector previously determined:
    highbin10 = posMx+o10dB;
    %Add something to remove values that are above the nyquist. bw3db is a
    %vector of NaNs to start with, so in this case we just skip the part of
    %saving those numbers and leave the NaNs.
    if highbin10 > length(f)%This means the whole bandwidth is not enough to reach the -10 dB mark
        %use the highest number.we want to add a large value, even if it's
        %not quite large enough
        high10dB = f(end);
    else
        high10dB = f(highbin10);
    end
    lowbin10 = posMx-e10dB;
    if lowbin10 <= 0 %This means the whole bandwidth is not enough to reach the -10 dB mark
        %use the lowest number. we want to add a large value, even if it's
        %not quite large enough
        low10dB = f(1);
    else
        low10dB = f(lowbin10);
    end
    bw10 = high10dB-low10dB;
    bw10db(c,:)= [low10dB, high10dB, bw10]; 
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Simone has code for calculating the SNR, which I"ll use for selecting only
    %the "best" clicks to look at. yNFilt and yFilt are the raw signal, before
    %fft or windowing or anything else, the same as mine. Also, yNFilt is
    %most likely the same values repeated over and over for each N clicks,
    %it's not N different noise samples (because of how the noise is
    %calculated), so just one value per set of clicks. 
    %calculate rms level of noise
    %assuming just one noise sample and one click, N = 1, so no need to do 1/N.
    yrms = sqrt(sum(yFilt{c,1}.*yFilt{c,1}));
    rmsSignal = 20*log10(yrms);
    %make them the same length.
    clickLength = size(yFilt{c,1},2);
    yNFiltCrop = yNFilt{c,1}(1:clickLength);
    yrmsNoise = sqrt(sum(yNFiltCrop.*yNFiltCrop));
    rmsNoise = 20*log10(yrmsNoise);
    
    snrTemp = rmsSignal - rmsNoise;
    snr(c) = snrTemp;

%     if snr > 15 %only looking at clicks with snr greater than 10
%         figure(1)
%         plot(linearSpec)
%         figure(2)
%         plot(f,specClickTf{c})
%     end
    
          %SBP Original code
% %       %calculate rms level of signal and noise
%         n = length(yNFilt);
%         yrmsNoise=zeros(size(yNFilt,1),1);
%         rmsNoise=zeros(size(yNFilt,1),1);
%         for i=1:size(yNFilt,1)
%             yrmsNoise(i) = sqrt(sum(yNFilt(i,:).*yNFilt(i,:))/n);
%             rmsNoise(i) = 20*log10(yrmsNoise(i));
%         end
% 
%         %calculate rms level of signal
%         %signal starts at tsFiltCountAll(101,:) and has duration(i)
%         yrms=zeros(size(yFilt,1),1);
%         rmsSignal=zeros(size(yFilt,1),1);
%         for i=1:size(yFilt,1)
%             n = round((dur(i)/1000)*fs);
%             if n>355
%                 n=355;
%             end
%             yrms(i) = sqrt(sum(yFilt(i,(101:101+n)).*yFilt(i,(101:101+n)))/n);
%             rmsSignal(i) = 20*log10(yrms(i));
%         end
% 
%         snr=rmsSignal-rmsNoise;
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%Little test for checking duration on 500 us window duration
%     if snr > 20 %only looking at clicks with snr greater than 10
%         figure(1)
%         plot(env5)
%         hold on
%         line([startbin startbin], [min(env5) max(env5)],'Color','r','LineWidth',1);
%         line([endbin endbin], [min(env5) max(env5)],'Color','r','LineWidth',1); 
%     end
%         




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate Q-value:
    %Q-3dB = the peakFr/-3dB bandwidth
    %Qrms = the centFr/rms bandwidth (Kyhn et al 2009)
    %Might as well get both, since other papers report both, and I have the
    %data for both. 
    QRMS_temp = peakFreq/bw3;
    QRMS(c) = QRMS_temp;
    Q3dB_temp = centroidFr/rmsBW;
    Q3dB(c) = Q3dB_temp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
dur95usec = dur95usec(clickInd,:);
dur95usecTails = dur95usecTails(clickInd,:);
durRMSus = durRMSus(clickInd,:);
bw3db = bw3db(clickInd,:);
bw10db = bw10db(clickInd,:);
bwRMS = bwRMS(clickInd,:);
QRMS = QRMS(clickInd,:);
Q3dB = Q3dB(clickInd,:);
% frames = frames{clickInd};
% yFilt = {yFilt{clickInd}};
yFilt = yFilt(clickInd);
yFiltBuff = yFiltBuff(clickInd);
specClickTf = specClickTf(clickInd);
specNoiseTf = specNoiseTf(clickInd);
peakFr = peakFr(clickInd,:);
centFr = centFr(clickInd,:);
snr = snr(clickInd,:);
%yNFilt = {yNFilt};
deltaEnv = deltaEnv(clickInd,:);
nDur = nDur(clickInd,:);
