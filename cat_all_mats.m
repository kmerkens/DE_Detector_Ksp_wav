%Code to combine the .mat files output from cat_click_times.m to get
%overall summary statistics and make plots of the detections from the
%entire deployment


close all
% clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Load data
% inDir = 'C:\Users\Karlina.Merkens\Documents\Kogia\DetectorOutput\Hawaii18';
% GraphDir = 'C:\Users\Karlina.Merkens\Documents\Kogia\DetectorOutput\Hawaii18\Final_histograms';
% inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\VJanik_Ksima_Wild\metadata\kogia';
% GraphDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\VJanik_Ksima_Wild\metadata\kogia\matlab_graphs';
% inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\DMann_Ksima_captive\metadata\kogia';
% GraphDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\DMann_Ksima_captive\metadata\kogia\matlab_graphs';
% inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\kogia\LowFreqTesting\metadata_100-119PeakFr\kogia';
% GraphDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\kogia\LowFreqTesting\metadata_100-119PeakFr\kogia\matlab_graphs';
% inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\metadata_170327_buzz\kogia';
% GraphDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\metadata_170327_buzz\kogia\matlab_graphs';
inDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\metadata\kogia';
GraphDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\NOAACRP_CNMI_Ksima_Wild\metadata\kogia\matlab_graphs';
% inpath = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\TGridley_Ksima_Wild\metadata\kogia';
% GraphDir = 'C:\Users\KMERKENS\Documents\Kogia\OtherRecordings\TGridley_Ksima_Wild\metadata\kogia\matlab_graphs';


% matList = dir(fullfile(inDir,'VJanik*.mat')); % Add wildcard to match the files you want to process.
% matList = dir(fullfile(inDir,'DMann*.mat')); % Add wildcard to match the files you want to process.
matList = dir(fullfile(inDir,'NOAA*.mat')); % Add wildcard to match the files you want to process.
% matList = dir(fullfile(inDir,'TGridley*.mat')); % Add wildcard to match the files you want to process.

%If you want to look at a subset of data based on the SNR, indicate that
%here, as well as the SNR threshold to use as a cutoff
subset = 0; %0 for no, 1 for yes
SNRthresh = 20; %If subset = 1 only clicks with SNR >= this value will be examined.

% fs = 375000; %Janik
% fs = 500000; %Mann
% fs = 384000; %CARB
% fs = 576000; %Gridley



allclickDnum = [];
alldurClickcon = [];
allnDurcon = [];
allndur95con = [];
allndur95Tailscon = [];
allbw3dbcon = [];
allbw10dbcon = [];
allbwRMScon = [];
allQRMScon = [];
allQ3dBcon = [];
allpeakFrcon = [];
allcentFrcon = [];
allsnrcon = [];
allppSignalcon = [];
allspecClickTfcon = [];
allspecNoiseTfcon = [];
allyFiltcon = [];
allmedianValues = [];
allmeanSpecClicks = [];
allmeanSpecNoises = [];
alliciEncs = [];

%%Concatenate data from all .mat files
for i1 = 1:length(matList)
    load(fullfile(inDir,matList(i1).name), 'clickDnum','durClickcon',...
        'nDurcon','ndur95con','ndur95Tailscon','bw3dbcon','bw10dbcon',...
        'bwRMScon','QRMScon','Q3dBcon','peakFrcon','centFrcon','snrcon','ppSignalcon','specClickTfcon',...
        'specNoiseTfcon','yFiltcon','medianValues','meanSpecClicks','meanSpecNoises',...
        'iciEncs','f','hdr')
    
    allclickDnum = [allclickDnum;clickDnum]; %save to one vector
    alldurClickcon = [alldurClickcon;durClickcon];
    allnDurcon = [allnDurcon; nDurcon];
    allndur95con = [allndur95con; ndur95con];
    allndur95Tailscon = [allndur95Tailscon; ndur95Tailscon];
    allbw3dbcon = [allbw3dbcon; bw3dbcon(:,3)];
    allbw10dbcon = [allbw10dbcon; bw10dbcon(:,3)];
    allbwRMScon = [allbwRMScon; bwRMScon];
    allQRMScon = [allQRMScon; QRMScon];
    allQ3dBcon = [allQ3dBcon; Q3dBcon];
    allpeakFrcon = [allpeakFrcon; peakFrcon];
    allcentFrcon = [allcentFrcon; centFrcon];
    allsnrcon = [allsnrcon; snrcon];
    allppSignalcon = [allppSignalcon; ppSignalcon];
    allspecClickTfcon = [allspecClickTfcon; specClickTfcon];
    allspecNoiseTfcon = [allspecNoiseTfcon; specNoiseTfcon];
    allyFiltcon = [allyFiltcon; yFiltcon];
    allmedianValues = [allmedianValues;medianValues];
    allmeanSpecClicks = [allmeanSpecClicks;meanSpecClicks];
    allmeanSpecNoises = [allmeanSpecNoises;meanSpecNoises];
    alliciEncs = [alliciEncs;iciEncs];
     
end

%convert click durations into us
allnDurcon = allnDurcon/(hdr.fs/1e6);

%%%%%
%Do you want to subset according to snr, to get the "best" clicks? Here's
%where you implement that subsetting. 
if subset == 1; %if yes
    trash = allsnrcon < SNRthresh; %Find the clicks that match this threshold. 
    
    allclickDnum(trash) = []; %save to one vector
    alldurClickcon(trash) = [];
    allnDurcon(trash) = [];
    allndur95con(trash) = [];
    allndur95Tailscon(trash) = [];
    allbw3dbcon(trash) = [];
    allbw10dbcon(trash) = [];
    allbwRMScon(trash) = [];
    allQRMScon(trash) = [];
    allQ3dBcon(trash) = [];
    allpeakFrcon(trash) = [];
    allcentFrcon(trash) = [];
    allsnrcon(trash) = [];
    allppSignalcon(trash) = [];
    allspecClickTfcon(trash) = [];
    allspecNoiseTfcon(trash) = [];
    allyFiltcon(trash) = [];
    %median values just emptied because they're not for this subset
    allmedianValues = [];
    allmeanSpecClicks(trash) = [];
    allmeanSpecNoises(trash) = [];
    %ICIs don't match the clicks, so keep them the same as before
    
end




%Ok, now you can proceed with the rest of the analysis
numclicks = size(allclickDnum,1);
strnclicks = num2str(numclicks);
numenc = size(allmeanSpecClicks,1);
strnumenc = num2str(numenc);

filedate = datestr(now, 'yymmdd');

%Add something to identify whether this applies to a subset of the data or
%the full data set. 
if subset == 1;
    fileappend = 'subset_';
else
    fileappend = '';
end

save([inDir,'\AllParamsConcat_',fileappend,filedate,'.mat'],...
'f','allclickDnum','alldurClickcon','allndur95con','allndur95Tailscon',...
'allbw3dbcon','allbw10dbcon','allbwRMScon','allQRMScon','allQ3dBcon',...
'allnDurcon','allpeakFrcon','allcentFrcon','allsnrcon',...
'allppSignalcon','allspecClickTfcon','allspecNoiseTfcon','allyFiltcon',...
'allmedianValues','allmeanSpecClicks','allmeanSpecNoises','alliciEncs')

%%%%
%Want to save everything in a .csv format? Then call this script:
saveDetecorOutput_170329(inDir,f,allndur95con,allbw3dbcon,...
    allbw10dbcon,allbwRMScon,allQRMScon,allQ3dBcon,...
    allpeakFrcon,allcentFrcon,allsnrcon,...
    allspecClickTfcon,allspecNoiseTfcon,...
    alliciEncs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate 25th, median and 75th percentiles for median values, to give some sense
%of width of the histogram spans. And also mean and stdev
Q1peakFr = prctile(allpeakFrcon,25);
Q1centFr = prctile(allcentFrcon,25);
Q1iciEncs = prctile(alliciEncs,25);
Q1durClick = prctile(alldurClickcon,25);
Q1ndur95 = prctile(allndur95con,25);
Q1ndur95Tails = prctile(allndur95Tailscon,25);
Q1bw3db = prctile(allbw3dbcon,25);
Q1bw10db = prctile(allbw10dbcon,25);
Q1bwRMS = prctile(allbwRMScon,25);
Q1QRMS = prctile(allQRMScon,25);
Q1Q3dB = prctile(allQ3dBcon,25);
Q1nDur = prctile(allnDurcon,25);
Q1ppSig = prctile(allppSignalcon,25);
Q1snr = prctile(allsnrcon,25);

medpeakFr = prctile(allpeakFrcon,50);
medcentFr = prctile(allcentFrcon,50);
mediciEncs = prctile(alliciEncs,50);
meddurClick = prctile(alldurClickcon,50);
medndur95 = prctile(allndur95con,50);
medndur95Tails = prctile(allndur95Tailscon,50);
medbw3db = prctile(allbw3dbcon,50);
medbw10db = prctile(allbw10dbcon,50);
medbwRMS = prctile(allbwRMScon,50);
medQRMS = prctile(allQRMScon,50);
medQ3dB = prctile(allQ3dBcon,50);
mednDur = prctile(allnDurcon,50);
medppSig = prctile(allppSignalcon,50);
medsnr = prctile(allsnrcon,50);

Q3peakFr = prctile(allpeakFrcon,75);
Q3centFr = prctile(allcentFrcon,75);
Q3iciEncs = prctile(alliciEncs,75);
Q3durClick = prctile(alldurClickcon,75);
Q3ndur95 = prctile(allndur95con,75);
Q3ndur95Tails = prctile(allndur95Tailscon,75);
Q3bw3db = prctile(allbw3dbcon,75);
Q3bw10db = prctile(allbw10dbcon,75);
Q3bwRMS = prctile(allbwRMScon,75);
Q3QRMS = prctile(allQRMScon,75);
Q3Q3dB = prctile(allQ3dBcon,75);
Q3nDur = prctile(allnDurcon,75);
Q3ppSig = prctile(allppSignalcon,75);
Q3snr = prctile(allsnrcon,75);


%Averages 
meandurClick = mean(alldurClickcon);
meannDur = mean(allnDurcon);
meanndur95 = mean(allndur95con);
meanndur95Tails = mean(allndur95Tailscon);
meanbw3db = mean(allbw3dbcon);
meanbw10db = mean(allbw10dbcon);
meanpeakFr = mean(allpeakFrcon);
meanppSignal = mean(allppSignalcon);
%meanspecClickTf = mean(allspecClickTfcon);
meaniciEncs = mean(alliciEncs);
meancentFr = mean(allcentFrcon);
meanbwRMS = mean(allbwRMScon);
meanQRMS = mean(allQRMScon);
meanQ3dB = mean(allQ3dBcon);
meansnr = mean(allsnrcon);


%Standard Deviation
stddurClick = std(alldurClickcon);
stdnDur = std(allnDurcon);
stdndur95 = std(allndur95con);
stdndur95Tails = std(allndur95Tailscon);
stdbw3db = std(allbw3dbcon);
stdbw10db = std(allbw10dbcon);
stdpeakFr = std(allpeakFrcon);
stdppSignal = std(allppSignalcon);
%stdspecClickTf = std(allspecClickTfcon);
stdiciEncs = std(alliciEncs);
stdcentFr = std(allcentFrcon);
stdbwRMS = std(allbwRMScon);
stdQRMS = std(allQRMScon);
stdQ3dB = std(allQ3dBcon);
stdsnr = std(allsnrcon);


%Standard Error - stdev / sqrt sample size
sedurClick = std(alldurClickcon)/sqrt(length(alldurClickcon));
sednDur = std(allnDurcon)/sqrt(length(allnDurcon));
sedndur95 = std(allndur95con)/sqrt(length(allndur95con));
sedndur95Tails = std(allndur95Tailscon)/sqrt(length(allndur95Tailscon));
sedbw3db = std(allbw3dbcon)/sqrt(length(allbw3dbcon));
sedbw10db = std(allbw10dbcon)/sqrt(length(allbw10dbcon));
sedpeakFr = std(allpeakFrcon)/sqrt(length(allpeakFrcon));
sedppSignal = std(allppSignalcon)/sqrt(length(allppSignalcon));
%sedspecClickTf = std(allspecClickTfcon)/sqrt(length(allspecClickTfcon));
sediciEncs = std(alliciEncs)/sqrt(length(alliciEncs));
sedcentFr = std(allcentFrcon)/sqrt(length(allcentFrcon));
sedbwRMS = std(allbwRMScon)/sqrt(length(allbwRMScon));
sedQRMS = std(allQRMScon)/sqrt(length(allQRMScon));
sedQ3dB = std(allQ3dBcon)/sqrt(length(allQ3dBcon));
sedsnr = std(allsnrcon)/sqrt(length(allsnrcon));



% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %Analysis of clicks short enough to only be singles. 
% % % [singleMeandurClick, singleMeddurClick] = singleclicksAnalysis(...
% % %     allclickDnum,alldurClickcon, allndur95con, allndur95Tailscon, ...
% % %     allbw3dbcon, allbw10dbcon, allbwRMScon, allnDurcon, allpeakFrcon, ...
% % %     allcentFrcon, allppSignalcon,allsnrcon,...
% % %     allspecClickTfcon,allspecNoiseTfcon, allyFiltcon,inDir,GraphDir,f);
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Make plots and save them
% 
% %Click duration - old method
% subplot(2,1,1)
% vec=0:10:400;
% hist(alldurClickcon,vec)
% xlim([0 400])
% line([meddurClick meddurClick], [0 200],'Color','r','LineWidth',1);
% line([Q1durClick Q1durClick], [0 200],'Color','r','LineWidth',1,'LineStyle',':');
% line([Q3durClick Q3durClick], [0 200],'Color','r','LineWidth',1,'LineStyle',':');
% text(0.5,0.9,['25th pctile = ',num2str(Q1durClick),' \musec'],'Unit','normalized','Color','r')
% text(0.5,0.8,['median = ',num2str(meddurClick),' \musec'],'Unit','normalized','Color','r')
% text(0.5,0.7,['75th pctile = ',num2str(Q3durClick),' \musec'],'Unit','normalized','Color','r')
% ylim([0 120])
% title(['Histgram of all click durations - old (\musec) (n = ',strnclicks,')']);
% ylabel('Counts')
% xlabel('Time (\musec)')
% 
% %nDur95 = Envelope duration 95%
% subplot(2,1,2)
% hist(allndur95con,vec)
% xlim([0 400])
% line([medndur95 medndur95], [0 200],'Color','r','LineWidth',1);
% line([Q1ndur95 Q1ndur95], [0 200],'Color','r','LineWidth',1,'LineStyle',':');
% line([Q3ndur95 Q3ndur95], [0 200],'Color','r','LineWidth',1,'LineStyle',':');
% text(0.1,0.9,['25th pctile = ',num2str(Q1ndur95),' \musec'],'Unit','normalized','Color','r')
% text(0.1,0.8,['median = ',num2str(medndur95),' \musec'],'Unit','normalized','Color','r')
% text(0.1,0.7,['75th pctile = ',num2str(Q3ndur95),' \musec'],'Unit','normalized','Color','r')
% ylim([0 70])
% title(['95% (\musec) (n = ',strnclicks,')']);
% %ylabel('Counts')
% xlabel('Time (\musec)')
% filename = fullfile(GraphDir,['Duration_',filedate]);
% saveas(gca, filename, 'tif')
% saveas(gca, filename, 'jpg')
% close


%Make plots and save them

% % % %Click duration
% % % %This version for Bahamas data
% % % %First, generate the subset of data with durClick shorter than 235 us
% % % shortdur = find(alldurClickcon > 235);
% % % shortset = alldurClickcon;
% % % shortset(shortdur) = [];
% % % figure(1)
% % % vec = (90:10:350);
% % % [n1, x1] = hist(alldurClickcon,vec);
% % % [n2, x2] = hist(shortset,vec);
% % % h1 = bar(x1, n1, 'hist');
% % % % set(h1, 'FaceColor','w','EdgeColor','k','LineWidth',2)
% % % set(h1,'FaceColor','k','EdgeColor','k','LineWidth',1);
% % % hold on
% % % h2 = bar(x2, n2,'hist');
% % % % set(h2,'FaceColor',[0.3,0.3,0.3],'EdgeColor',[0.3,0.3,0.3]);
% % % set(h2, 'FaceColor',[0.5,0.5,0.5],'EdgeColor',[0.5,0.5,0.5],'LineWidth',1)
% % % % line([meddurClick meddurClick], [0 50],'Color','k','LineWidth',2,'LineStyle','--'); %median of all clicks 184 us
% % % % line([meandurClick meandurClick], [0 50],'Color','k','LineWidth',2,'LineStyle','--'); %mean of all clicks,
% % % % line([singleMeddurClick singleMeddurClick], [0 50],'Color','k','LineWidth',2); %median of the single clicks,162 us
% % % % line([singleMeandurClick singleMeandurClick], [0 50],'Color','k','LineWidth',2); %mean of the single clicks
% % % line([meddurClick meddurClick], [0 50],'Color','k','LineWidth',2); %median of all clicks 184 us
% % % line([meandurClick meandurClick], [0 50],'Color','k','LineWidth',2); %mean of all clicks,
% % % line([singleMeddurClick singleMeddurClick], [0 50],'Color','k','LineWidth',2,'LineStyle','--'); %median of the single clicks,162 us
% % % line([singleMeandurClick singleMeandurClick], [0 50],'Color','k','LineWidth',2,'LineStyle','--'); %mean of the single clicks
% % % % text(0.01,0.9,['Median of all clicks = ',num2str(meddurClick),' \musec'],'Unit','normalized','Color','k')
% % % % text(0.01,0.85,['Median of clicks shorter than 235 \musec = ',num2str(singleMeddurClick),' \musec'],'Unit','normalized','Color','k')
% % % xlim([80 360])
% % % ylim([0 50])
% % % % title(['Histgram of all click durations (\musec) (n = ',strnclicks,')']);
% % % ylabel('Counts','FontSize',18)
% % % xlabel('Time (\musec)','FontSize',18)
% % % set(gca,'YTick',[0:10:45],'FontSize',18);
% % % %set shape/size
% % % set(gcf,'PaperUnits', 'inches','PaperPosition', [0 0 6 5]);
% % % filename = fullfile(GraphDir,['SingleClick_MeanSpecGrand_forPub',filedate]);
% % % print('-dtiff',filename)
% % % %print('-djpeg', filename)
% % % print('-dpng', filename,'-r600')
% % % filename = fullfile(GraphDir,['Click_duration_pretty_8_',filedate]);
% % % saveas(gca, filename, 'tif')
% % % saveas(gca, filename, 'jpg')
% % % close(figure(1));

%Click duration
figure(1)
hist(alldurClickcon,50)
line([meddurClick meddurClick], [0 60],'Color','r','LineWidth',3);
line([Q1durClick Q1durClick], [0 60],'Color','r','LineWidth',2,'LineStyle',':');
line([Q3durClick Q3durClick], [0 60],'Color','r','LineWidth',2,'LineStyle',':');
text(0.5,0.9,['25th pctile = ',num2str(Q1durClick,3),' \musec'],'Unit','normalized','Color','r')
text(0.5,0.85,['median = ',num2str(meddurClick,3),' \musec'],'Unit','normalized','Color','r')
text(0.5,0.8,['75th pctile = ',num2str(Q3durClick,3),' \musec'],'Unit','normalized','Color','r')
% ylim([0 5])
title(['Histgram of click durations (\musec) (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('Time (\musec)')
filename = fullfile(GraphDir,['Click_duration_',fileappend,filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
close(figure(1));



% %nDur = Envelope duration
% figure(2)
% hist(allnDurcon,25)
% line([mednDur mednDur], [0 700],'Color','r','LineWidth',3);
% line([Q1nDur Q1nDur], [0 700],'Color','r','LineWidth',2,'LineStyle',':');
% line([Q3nDur Q3nDur], [0 700],'Color','r','LineWidth',2,'LineStyle',':');
% text(0.6,0.9,['25th pctile = ',num2str(Q1nDur,3),' \musec'],'Unit','normalized','Color','r')
% text(0.6,0.85,['median = ',num2str(mednDur,3),' \musec'],'Unit','normalized','Color','r')
% text(0.6,0.8,['75th pctile = ',num2str(Q3nDur,3),' \musec'],'Unit','normalized','Color','r')
% %ylim([0 1050])
% title(['Histgram of all click envelope durations (\musec) (n = ',strnclicks,')']);
% ylabel('Counts')
% xlabel('Time (\musec)')
% filename = fullfile(GraphDir,['Envelope_duration_',filedate]);
% saveas(gca, filename, 'tif')
% saveas(gca, filename, 'jpg')
% close(figure(2));

% 
%PeakFrequency
figure(3)
vec=0:1:140;
% hist(allpeakFrcon,40)
hist(allpeakFrcon,vec)
% line([medpeakFr medpeakFr], [0 1600],'Color','r','LineWidth',3);
% line([Q1peakFr Q1peakFr], [0 1600],'Color','r','LineWidth',2,'LineStyle',':');
% line([Q3peakFr Q3peakFr], [0 1600],'Color','r','LineWidth',2,'LineStyle',':');
% text(0.1,0.9,['25th pctile = ',num2str(Q1peakFr,3),' kHz'],'Unit','normalized','Color','r')
% text(0.1,0.85,['median = ',num2str(medpeakFr,3),' kHz'],'Unit','normalized','Color','r')
% text(0.1,0.8,['75th pctile = ',num2str(Q3peakFr,3),' kHz'],'Unit','normalized','Color','r')
%ylim([0 600])
%xlim([105 155])
% title(['Histgram of all peak frequencies (kHz) (n = ',strnclicks,')']);
% ylabel('Counts')
% xlabel('Frequency (kHz)')
%%%Making pretty for publication
xlim([100 135])
set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
set(findobj(gca,'Type','text'),'FontName','Arial','FontSize',12)
set(gca,'FontSize',14,'FontName','Arial')
ylabel('Counts')
xlabel('Frequency (kHz)')
title(['Histogram of peak frequencies (kHz)']);

set(gcf,'PaperUnits', 'inches','PaperPosition', [0 0 7.5 10]);
filename = fullfile(GraphDir,['Peak_Frequency_',fileappend,filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
print('-dpng', filename,'-r600')
close(figure(3));

%Centroid Frequency
figure(3)
vec=0:1:140;
hist(allcentFrcon,vec)
%%%Making pretty for publication
xlim([100 135])
set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
set(findobj(gca,'Type','text'),'FontName','Arial','FontSize',12)
set(gca,'FontSize',14,'FontName','Arial')
ylabel('Counts')
xlabel('Frequency (kHz)')
title(['Histogram of centroid frequencies (kHz)']);

set(gcf,'PaperUnits', 'inches','PaperPosition', [0 0 7.5 10]);
filename = fullfile(GraphDir,['Cent_Frequency_',fileappend,filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
print('-dpng', filename,'-r600')
close(figure(3));


% 
% 
%ppSignal
figure(4)
hist(allppSignalcon,50)
line([medppSig medppSig], [0 80],'Color','r','LineWidth',3);
line([Q1ppSig Q1ppSig], [0 80],'Color','r','LineWidth',2,'LineStyle',':');
line([Q3ppSig Q3ppSig], [0 80],'Color','r','LineWidth',2,'LineStyle',':');
text(0.6,0.9,['25th pctile = ',num2str(Q1ppSig,3),' dB'],'Unit','normalized','Color','r')
text(0.6,0.85,['median = ',num2str(medppSig,3),' dB'],'Unit','normalized','Color','r')
text(0.6,0.8,['75th pctile = ',num2str(Q3ppSig,3),' dB'],'Unit','normalized','Color','r')
%ylim([0 5])
%xlim([105 155])
title(['Histgram of all peak-peak amplitudes (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('p-p amplitude (dB)')
filename = fullfile(GraphDir,['P-P_Amplitude_',fileappend,filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
close(figure(4));


%ici per encounter
%Remove ICIs less than 2 and greater than 400
remove = find(alliciEncs<2 | alliciEncs>390); %for usual clicks
% remove = find(alliciEncs<2 | alliciEncs>70); %To remove any ICIs that are really missed clicks for buzzes
alliciEncs(remove) = [];
numicis = size(alliciEncs,1);
strnicis = num2str(numicis);
figure(5)
 xbins = 10:5:405;
xbins = 0:2:100;
hist(alliciEncs,xbins)
% hist(alliciEncs,50)
% line([mediciEncs mediciEncs], [0 900],'Color','r','LineWidth',3);
% line([Q1iciEncs Q1iciEncs], [0 900],'Color','r','LineWidth',2,'LineStyle',':');
% line([Q3iciEncs Q3iciEncs], [0 900],'Color','r','LineWidth',2,'LineStyle',':');
% text(0.5,0.9,['25th pctile = ',num2str(Q1iciEncs,3),' msec'],'Unit','normalized','Color','r')
% text(0.5,0.85,['median = ',num2str(mediciEncs,3),' msec'],'Unit','normalized','Color','r')
% text(0.5,0.8,['75th pctile = ',num2str(Q3iciEncs,3),' msec'],'Unit','normalized','Color','r')
% ylim([0 40])
% xlim([20 450]) 
xlim([0 100]) %change to 100 for buzzes since I cut off everything above 1100 ms
title(['Histgram of all inter-click-intervals (n = ',strnicis,')']);
ylabel('Counts')
xlabel('Time (msec)')
set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
filename = fullfile(GraphDir,['ICI_',fileappend,filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
close(figure(5));




%-3 dB Bandwidth
figure(5)
xbins = 0:1:35;
hist(allbw3dbcon,xbins)
xlim([0 35])
title(['Histgram of all -3dB Bandwidths (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('-3dB Bandwidth (kHz)')
set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
filename = fullfile(GraphDir,['3dBBW_',fileappend,filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
% saveas(gca, filename, 'fig')
close(figure(5));

%-10 dB Bandwidth
figure(5)
% xbins = 0:2:130;
xbins = 0:1:35;
hist(allbw10dbcon,xbins)
xlim([0 35])
title(['Histgram of all -10dB Bandwidths (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('-10dB Bandwidth (kHz)')
set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
filename = fullfile(GraphDir,['10dBBW_',fileappend,filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
% saveas(gca, filename, 'fig')
close(figure(5));

%RMS Bandwidth
figure(6)
xbins = 0:1:45;
hist(allbwRMScon,xbins)
xlim([0 45])
% ylim([0 100])
title(['Histgram of all RMS Bandwidths (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('RMS Bandwidth (kHz)')
set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
filename = fullfile(GraphDir,['RMSBW_',fileappend,filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
% saveas(gca, filename, 'fig')
close(figure(6));

%SNR
figure(7)
xbins = 0:1:30;
hist(allsnrcon,xbins)
xlim([0 30])
title(['Histgram of all SNR (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('SNR')
set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
filename = fullfile(GraphDir,['SNR_',fileappend,filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
% saveas(gca, filename, 'fig')
close(figure(7));


% % %Next -plot the median values
% % numMeds = size(allmedianValues,1);
% % strnMeds = num2str(numMeds);
% % titlecell = {'Peak Frequency', 'Inter-click-interval','Click Duration',...
% %     'Peak-to-peak amplitude';'(kHz)','(ms)','(\musec)', '(dB)';...
% %     'Frequency', 'Time','Time','Amplitude'};
% % figure(6)
% % for p = 1:4
% %     subplot(2,2,p)
% %     hist(allmedianValues(:,p),30);
% %     meanMed = nanmean(allmedianValues(:,p));
% %     text(0.1,0.9,['mean = ',num2str(meanMed,3),' ',titlecell{2,p}],'Unit','normalized','Color','r')
% %     title(titlecell{1,p});
% %     xlabel([titlecell{3,p},' ',titlecell{2,p}])
% %     ylabel('Counts')
% % end
% % text(-1,2.55,['Histgrams of median parameters per Encounter (n = ',strnMeds,')'],'Unit','normalized')
% % filename = fullfile(GraphDir,['MediansPerEncounter',filedate]);
% % saveas(gca, filename, 'tif')
% % saveas(gca, filename, 'jpg')
% % close(figure(6));


%Then plot the mean spectra - could overlay with different colors
numSpecs = size(allmeanSpecClicks,1);
strnSpecs = num2str(numSpecs);
%Make frequency axis
numFreqBins = (size(allmeanSpecClicks,2))-1;
%steps = (fs/2)/(numFreqBins-1);
%freqs = (0:steps:(fs/2))/1000;
freqs = f;

% figure(7)
% for p = 1:numSpecs
%     %plot(freqs,allmeanSpecClicks(p,2:end), 'Color',rand(1,3));
%     plot(freqs,allmeanSpecClicks(p,2:end), 'Color','k','LineWidth',1);
%     hold on
%     plot(freqs,allmeanSpecNoises(p,2:end), 'Color',[0.7, 0.7, 0.7], 'LineWidth',1);
% end
% legend('mean click','mean noise','Location','southwest')
% title(['Mean Spectrum of Each Encounter (n = ',strnumenc,')']);
% xlabel('Frequency (kHz)')
% ylabel('Amplitude (dB)')
% filename = fullfile(GraphDir,['MeanSpecPerEncounter',filedate]);
% saveas(gca, filename, 'tif')
% saveas(gca, filename, 'jpg')
% close(figure(7));
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Then combine all spectra and calculate one mean, plot with the mean noise.
%First, calculate the mean. start by looping through all the entries in the
%cell array and saving those to a non-cell matrix
numspecClick = size(allspecClickTfcon,1);
numspecNoise = size(allspecNoiseTfcon,1);
concatspecs = [];
concatspecsN = [];
for s = 1:numspecClick
        concatspecs = [concatspecs;allspecClickTfcon{s,1}'];
end
grandmeanSpec = mean(concatspecs);


for s = 1:numspecNoise
        concatspecsN = [concatspecsN;allspecNoiseTfcon{s,1}'];
end
%Get unique noise spectra (because some/many are repeats)
concatspecsNU = unique(concatspecsN,'rows');
grandmeanSpecN = mean(concatspecsNU);

%Make frequency axis
numFreqBins = (size(allmeanSpecClicks,2))-1;
freqs = f;

    
figure(8)
plot(freqs,grandmeanSpec(1,1:end), 'Color','k','LineWidth',1);
hold on
plot(freqs,grandmeanSpecN(1,1:end), 'Color',[0.7, 0.7, 0.7], 'LineWidth',1);
% plot(freqs,allmeanSpecNoises(1,2:end), 'Color',[0.7, 0.7, 0.7], 'LineWidth',1);
xlim([0 175]);
legend('mean click','mean noise','Location','southwest')
title(['Mean Spectrum of all Clicks (n = ',num2str(numspecClick),')']);
xlabel('Frequency (kHz)')
ylabel('Amplitude (dB)')
filename = fullfile(GraphDir,['MeanSpecGrand',fileappend,filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
close(figure(8));


%Make it pretty for publications
%Normalize to max value
prenorm = grandmeanSpec(1,1:end);
maxfreq = max(prenorm);
normdbs = prenorm-maxfreq;
prenormnoise = grandmeanSpecN(1,1:end);
normnoisedbs = prenormnoise-maxfreq;
%Truncate at 20kHz and 170kHz 
toolowfreq = find(freqs < 20);
normdbs(toolowfreq) = [];
normnoisedbs(toolowfreq) = [];
plotfreqs = freqs;
plotfreqs(toolowfreq) = [];
toohighfreq = find(plotfreqs > 180); %180 for all since BPF at 170 now
normdbs(toohighfreq) = [];
normnoisedbs(toohighfreq) = [];
plotfreqs(toohighfreq) = [];

figure(9)
plot(plotfreqs,normdbs, 'Color','k','LineWidth',2);
hold on
plot(plotfreqs,normnoisedbs, 'Color',[0.7, 0.7, 0.7],...
    'LineWidth',2);
xlim([0 175]); %original 200. Up to 230 for D mann
set(gca,'XTick',[0,40,80,120,160,200],'FontSize',18);
xlabel('Frequency (kHz)','FontSize',18);
ylim([-68 3]); %oroginal 42, up to 45 for D Mann, 68 for TGridley
% set(gca,'YTick',[-30,-15,0],'FontSize',18); %For most
set(gca,'YTick',[-60:10:0],'FontSize',18); %For TGridley
ylabel('Relative Spectral Level (dB)','FontSize',18);
%set shape/size
set(gcf,'PaperUnits', 'inches','PaperPosition', [0 0 6 6]);

filename = fullfile(GraphDir,['MeanSpecGrand_forPub',fileappend,filedate]);
print('-dtiff',filename)
%print('-djpeg', filename)
print('-dpng', filename,'-r600')
%print('-deps', filename,'-r200')
% saveas(gca, filename, 'fig')
close(figure(9));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
% %make a time series!
% figure(9)
% hist(allclickDnum(:,1),77) %77 days between start and end of effort, to get clicks per day.
% ylim([0 1300])
% % xlim([0 0.5])
% %Add grey boxes for no refording effort. Lines at july 27 and october 11 
% dategapblocksX = [735781 735781 735807 735807 735885 735885 735903 735903];%values for bottoms and tops of squares to add [b,t,t,b]
% dategapblocksY = [0 1300 1300 0 0 1300 1300 0];
% hold on
% grey = [0.8,0.8,0.8];
% fill(dategapblocksX,dategapblocksY, grey);
% ylim([0 1300])
% current_day = datenum('Jul 15 2014'); %the date the plot starts
% [Y, M, D, H, MN, S] = datevec(current_day);
% current_day = addtodate(current_day, -D + 1,'day');
% last_day = datenum('Oct 31 2014');
% xlim([current_day last_day]);
% xtick = [current_day];
% xstep_length = 1; %one tick per month
% xstep_unit = 'month';
% while (last_day > current_day)
%     xtick = [xtick addtodate(current_day, xstep_length, xstep_unit)];
%     current_day = addtodate(current_day, xstep_length, xstep_unit);
% end
% ts = gca;
% set(ts,'XTick', xtick,'XTickLabel', datestr(xtick,'mmm yy')); 
% %,'FontSize',14
% set(gca, 'Ticklength', [0 0])
% set(gcf, 'renderer', 'zbuffer'); %removes the exponent from the date  
% title(['Time Series of All Encounters (n = ',strnumenc,')']);
% ylabel('Counts of Clicks per Day')
% xlabel('Date') %Need something here to make this axis pretty with months/dates
% filename = fullfile(GraphDir,['TS_',filedate]);
% saveas(gca, filename, 'tif')
% saveas(gca, filename, 'jpg')
% close(figure(9));



%summary stats not used.
% %Skewness
% skedurClick = skewness(alldurClickcon);
% skenDur = skewness(allnDurcon);
% skendur95 = skewness(allndur95con);
% skendur95Tails = skewness(allndur95Tailscon);
% skebw3db = skewness(allbw3dbcon);
% skebw10db = skewness(allbw10dbcon);
% skepeakFr = skewness(allpeakFrcon);
% skeppSignal = skewness(allppSignalcon);
% %skespecClickTf = skewness(allspecClickTfcon);
% skeiciEncs = skewness(alliciEncs);



%Make one table, save it as .mat and .xls 
SummaryStats = [Q1durClick Q1ndur95 Q1ndur95Tails Q1bw3db Q1bw10db Q1bwRMS Q1QRMS Q1Q3dB Q1peakFr Q1centFr Q1ppSig Q1iciEncs Q1snr;...
    meandurClick meanndur95 meanndur95Tails meanbw3db meanbw10db meanbwRMS meanQRMS meanQ3dB meanpeakFr meancentFr meanppSignal meaniciEncs meansnr;...
    meddurClick medndur95 medndur95Tails medbw3db medbw10db medbwRMS medQRMS medQ3dB medpeakFr medcentFr medppSig mediciEncs medsnr;...
    Q3durClick Q3ndur95 Q3ndur95Tails Q3bw3db Q3bw10db Q3bwRMS Q3QRMS Q3Q3dB Q3peakFr Q3centFr Q3ppSig Q3iciEncs Q3snr;...
    stddurClick stdndur95 stdndur95Tails stdbw3db stdbw10db stdbwRMS stdQRMS stdQ3dB stdpeakFr stdcentFr stdppSignal stdiciEncs stdsnr;...
    sedurClick sedndur95 sedndur95Tails sedbw3db sedbw10db sedbwRMS sedQRMS sedQ3dB sedpeakFr sedcentFr sedppSignal sediciEncs sedsnr];
    %skedurClick skendur95 skendur95Tails skebw3db skebw10db skepeakFr skeppSignal skeiciEncs];


filename = fullfile(inDir,['SummaryStats_',fileappend,filedate]);
save(filename, 'SummaryStats')
xlswrite([filename,'.xls'],SummaryStats)
%Write csv with headers
headers = {'durClick','ndur95','ndur95Tails','bw3db','bw10db',...
    'bwRMS','QRMS','Q3dB','peakFr','centFr','ppSig','iciEncs','snr'};
csvwrite_with_headers([filename,'.csv'],SummaryStats,headers)
