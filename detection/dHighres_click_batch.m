function dHighres_click_batch(fullFiles,fullLabels,inDisk,p,viewPath,...
    tfFullFile,encounterTimes,guideDetector,GraphDir,DASBR)

N = length(fullFiles);
previousFs = 0; % make sure we build filters on first pass

for idx1=1:N; % for each data file
    
    % figure out which files are needed, where to find them.
    [hdr,channel,labelFile]...
        = dInput_HR_files(fullFiles{idx1},fullLabels{idx1},viewPath,p);
    
    if isempty(hdr.fs)
        continue % skip if you couldn't read a header
    elseif hdr.fs ~= previousFs
        % otherwise, if this is the first time through, build your filters,
        % only need to do this once though, so if you already have this
        % info, this step is skipped
        [previousFs,fftSize,fftWindow,binWidth_Hz,freq_kHz,...
            wideBandFilter,specRange] = dBuild_filters(p,hdr.fs);
        
        % Determine the frequencies for which we need the transfer function
        xfr_f = (specRange(1)-1)*binWidth_Hz:binWidth_Hz:(specRange(end)-1)*binWidth_Hz;
        if ~isempty(tfFullFile)
            [xfr_f, xfrOffset] = dtf_map(tfFullFile, xfr_f);
        else
            % if you didn't provide a tf function, then just create a
            % vector of zeros of the right size.
            xfrOffset = zeros(size(xfr_f));
        end
        xfrOffset = xfrOffset';
    end
    
    if exist([inDisk 'metadata\' labelFile],'file')
        % Read in the .c file produced by the short term detector.
        [starts,stops,~] = ioReadLabelFile([inDisk 'metadata\' labelFile]);
    end
    % Open xwav file
    fid = ioOpenViewpath(fullFiles{idx1}, viewPath, 'r');
    
    % Look for clicks, hand back parameters of retained clicks
    [clickTimes,ppSignalVec,durClickVec,dur95Vec,dur95TailsVec,durRMSusVec,bw3dbVec,...
        bw10dbVec,bwRMSVec,QRMSVec,Q3dBVec,~,yFiltVec,specClickTfVec,specNoiseTfVec,...
        peakFrVec,centFrVec,snrVec,yFiltBuffVec,f,deltaEnvVec,nDurVec]...
        = dProcess_HR_starts(fid, wideBandFilter,starts,stops,channel,...
        xfrOffset,specRange,p,hdr,fullFiles{idx1},fftWindow,fullLabels{idx1});
    
    % Done with that file
    fclose(fid);
    fclose all;
    %fprintf('done with %s\n', fullFiles{idx1});
    %toc
    
    % Run post processing to remove rogue loner clicks, prior to writing
    % the remaining output files.
    %Also, if DASPR remove any clicks from the first 0.06 seconds
    delFlag = clickInlinePProc(fullLabels{idx1},clickTimes,p,encounterTimes,...
        guideDetector,hdr,DASBR);
    delIdx = find(delFlag==1);
    
    % save a mat file now, rather than recalculating later
    clickTimes = clickTimes(delIdx,:);
    ppSignal = ppSignalVec(delIdx,:);
    durClick = durClickVec(delIdx,:);
    dur95 = dur95Vec(delIdx,:);
    dur95Tails = dur95TailsVec(delIdx,:);
    durRMSus = durRMSusVec(delIdx,:);
    bw3db = bw3dbVec(delIdx,:);
    bw10db = bw10dbVec(delIdx,:);
    bwRMS = bwRMSVec(delIdx,:);
    QRMS = QRMSVec(delIdx,:);
    Q3dB = Q3dBVec(delIdx,:);
    % yNFilt = yNFiltVec;
    specClickTf = specClickTfVec(delIdx,:);
    specNoiseTf = specNoiseTfVec(delIdx,:);
    peakFr = peakFrVec(delIdx,:);
    centFr = centFrVec(delIdx,:);
    snr = snrVec(delIdx,:);
    deltaEnv = deltaEnvVec(delIdx,:);
    nDur = nDurVec(delIdx,:);
    
    if ~isempty(delIdx)
        yFilt = yFiltVec(delIdx);
        yFiltBuff = yFiltBuffVec(delIdx);
    else
        yFilt = {};
        yFiltBuff = {};
    end
    
    save(strcat(fullLabels{idx1}(1:end-2),'.mat'),'clickTimes','ppSignal',...
        'durClick','dur95','dur95Tails','durRMSus',...
        'bw3db','bw10db','bwRMS','QRMS','Q3dB','f','hdr','nDur','deltaEnv','yFilt','specClickTf', ...
        'specNoiseTf','peakFr','centFr','snr','-mat','yFiltBuff');%
    
%     % Make plots of each encounter's click parameters, now that all extra
%     % clicks have been removed. Save those to a directory.
%     if guideDetector == 1
%         plotClickEncounters_150310(encounterTimes,clickTimes,ppSignal,...
%             durClick,dur95,bw3db,bw10db,specClickTf,specNoiseTf,peakFr,nDur,yFilt,hdr,GraphDir,f);
%     end
end
    

    
    
    
    
