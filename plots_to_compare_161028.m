%Code to open more than one "grand mean" data set from the click detector,
%and make comparison plots of the data. 
%To be run AFTER cat_all_mats.m
inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\kogia\ComparePlot\Buzz_Usual';
GraphDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\kogia\ComparePlot\Buzz_Usual';
matList = dir(fullfile(inDir,'AllParamsConcat*.mat')); % Add wildcard to match the files you want to process.

% fs = 375000; %Janik
% fs = 500000; %Mann
fs = 384000; %CARB

for i1 = 1:length(matList)
    load(fullfile(inDir,matList(i1).name),...
    'allclickDnum','alldurClickcon','allbw3dbcon','allbw10dbcon','allnDurcon','allpeakFrcon',...
    'allppSignalcon','allspecClickTfcon','allspecNoiseTfcon','allyFiltcon',...
    'allmedianValues','allmeanSpecClicks','allmeanSpecNoises','alliciEncs');

    %Save all the parameters in a cell array for later access
    clickDnum(i1,:) = {allclickDnum};
    durClick(i1,:) = {alldurClickcon};
    bw3db(i1,:) = {allbw3dbcon};
    bw10db(i1,:) = {allbw10dbcon};
    nDur(i1,:) = {allnDurcon};
    peakFr(i1,:) = {allpeakFrcon};
    ppSignal(i1,:) = {allppSignalcon};
    specClickTf(i1,:) = {allspecClickTfcon};
    specNoiseTf(i1,:) = {allspecNoiseTfcon};
    iciEncs(i1,:) = {alliciEncs};
end

   
    numenc = size(allmeanSpecClicks,1);
    strnumenc = num2str(numenc);

    filedate = datestr(now, 'yymmdd');
    
    %Make plots
    
    %Click duration
    figure(1)
    for p = 1:(size(durClick,1))
        hist(durClick{p},20,'Color',rand(1:3))
        numclicks = size(durClick{p},1);
        strnclicks = num2str(numclicks);
        hold on
    end

    title(['Histgram of all click durations (\musec)']);
    %(n = ',strnclicks,')']);
    ylabel('Counts')
    xlabel('Time (\musec)')
    filename = fullfile(GraphDir,['Click_duration_',filedate]);
    saveas(gca, filename, 'tif')
    saveas(gca, filename, 'jpg')
    close(figure(1));


    %PeakFrequency
    figure(2)
    hist(allpeakFrcon,25)
    line([medpeakFr medpeakFr], [0 1600],'Color','r','LineWidth',3);
    line([Q1peakFr Q1peakFr], [0 1600],'Color','r','LineWidth',2,'LineStyle',':');
    line([Q3peakFr Q3peakFr], [0 1600],'Color','r','LineWidth',2,'LineStyle',':');
    text(0.1,0.9,['25th pctile = ',num2str(Q1peakFr,3),' kHz'],'Unit','normalized','Color','r')
    text(0.1,0.85,['median = ',num2str(medpeakFr,3),' kHz'],'Unit','normalized','Color','r')
    text(0.1,0.8,['75th pctile = ',num2str(Q3peakFr,3),' kHz'],'Unit','normalized','Color','r')
    %ylim([0 600])
    %xlim([105 155])
    title(['Histgram of all peak frequencies (kHz) (n = ',strnclicks,')']);
    ylabel('Counts')
    xlabel('Frequency (kHz)')
    filename = fullfile(GraphDir,['Peak_Frequency_',filedate]);
    saveas(gca, filename, 'tif')
    saveas(gca, filename, 'jpg')
    %close(figure(2));


    %ici per encounter
    numicis = size(alliciEncs,1);
    strnicis = num2str(numicis);
    figure(3)
    hist(alliciEncs,50)
    line([mediciEncs mediciEncs], [0 900],'Color','r','LineWidth',3);
    line([Q1iciEncs Q1iciEncs], [0 900],'Color','r','LineWidth',2,'LineStyle',':');
    line([Q3iciEncs Q3iciEncs], [0 900],'Color','r','LineWidth',2,'LineStyle',':');
    text(0.5,0.9,['25th pctile = ',num2str(Q1iciEncs,3),' msec'],'Unit','normalized','Color','r')
    text(0.5,0.85,['median = ',num2str(mediciEncs,3),' msec'],'Unit','normalized','Color','r')
    text(0.5,0.8,['75th pctile = ',num2str(Q3iciEncs,3),' msec'],'Unit','normalized','Color','r')
    %ylim([0 375])
    % xlim([0 0.5])
    title(['Histgram of all inter-click-intervals (n = ',strnicis,')']);
    ylabel('Counts')
    xlabel('Time (msec)')
    filename = fullfile(GraphDir,['ICI_',filedate]);
    saveas(gca, filename, 'tif')
    saveas(gca, filename, 'jpg')
    %close(figure(3));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Then combine all spectra and calculate one mean, plot with the mean noise
    %from the first file.
    %First, calculate the mean. start by looping through all the entries in the
    %cell array and saving those to a non-cell matrix
    numspecClick = size(allspecClickTfcon,1);
    concatspecs = [];
    for s = 1:numspecClick
            concatspecs = [concatspecs;allspecClickTfcon{s,1}'];
    end
    grandmeanSpec = mean(concatspecs);

    %Make frequency axis
    numFreqBins = (size(allmeanSpecClicks,2))-1;
    freqs = f;

    figure(4)
    plot(freqs,grandmeanSpec(1,1:end), 'Color','k','LineWidth',1);
    hold on
    plot(freqs,allmeanSpecNoises(1,2:end), 'Color',[0.7, 0.7, 0.7], 'LineWidth',1);
    xlim([0 200]);
    legend('mean click','mean noise','Location','southwest')
    title(['Mean Spectrum of all Clicks (n = ',num2str(numspecClick),')']);
    xlabel('Frequency (kHz)')
    ylabel('Amplitude (dB)')
    filename = fullfile(GraphDir,['MeanSpecGrand',filedate]);
    saveas(gca, filename, 'tif')
    saveas(gca, filename, 'jpg')
    %close(figure(4));



