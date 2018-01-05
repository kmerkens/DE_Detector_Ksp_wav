%Testing the characteristics of K sima clicks from The Bahamas with
%duration less than 235 us (which are likely to just be single clicks, no
%echoes). ICIs are not associated with individual clicks, so that parameter
%is not considered here. 
 
function [meandurClick, meddurClick] = singleclicksAnalysis(...
    allclickDnum,alldurClickcon, allndur95con, allndur95Tailscon, ...
    allbw3dbcon, allbw10dbcon, allbwRMScon, allnDurcon, allpeakFrcon,...
    allcentFrcon, allppSignalcon,allsnrcon,...
    allspecClickTfcon,allspecNoiseTfcon, allyFiltcon,inDir,GraphDir,f);

%First identify and cut out the clicks that have duration longer than
%235us. 
echoclicks = alldurClickcon > 235; %Changed from > to <= to test second mode

singleclickDnum = [allclickDnum]; %save to one vector
singledurClickcon = [alldurClickcon];
singlenDurcon = [allnDurcon];
singlendur95con = [allndur95con];
singlendur95Tailscon = [allndur95Tailscon];
singlebw3dbcon = [allbw3dbcon];
singlebw10dbcon = [allbw10dbcon];
singlebwRMScon = [allbwRMScon];
singlepeakFrcon = [allpeakFrcon];
singlecentFrcon = [allcentFrcon];
singleppSignalcon = [allppSignalcon];
singlespecClickTfcon = [allspecClickTfcon];
singlespecNoiseTfcon = [allspecNoiseTfcon];
singleyFiltcon = [allyFiltcon];
singlesnrcon = [allsnrcon];


singleclickDnum(echoclicks,:) = []; %save to one vector
singledurClickcon(echoclicks) = [];
singlenDurcon(echoclicks) = [];
singlendur95con(echoclicks) = [];
singlendur95Tailscon(echoclicks) = [];
singlebw3dbcon(echoclicks) = [];
singlebw10dbcon(echoclicks) = [];
singlebwRMScon(echoclicks) = [];
singlepeakFrcon(echoclicks) = [];
singlecentFrcon(echoclicks) = [];
singleppSignalcon(echoclicks) = [];
singlespecClickTfcon(echoclicks) = [];
singlespecNoiseTfcon(echoclicks) = [];
singleyFiltcon(echoclicks) = [];
singlesnrcon(echoclicks) = [];


numclicks = size(singleclickDnum,1);
strnclicks = num2str(numclicks);
filedate = datestr(now, 'yymmdd');



%Then calculate the means, medians and quantiles for this set of clicks. 
Q1peakFr = prctile(singlepeakFrcon,25);
Q1centFr = prctile(singlecentFrcon,25);
Q1durClick = prctile(singledurClickcon,25);
Q1ndur95 = prctile(singlendur95con,25);
Q1ndur95Tails = prctile(singlendur95Tailscon,25);
Q1bw3db = prctile(singlebw3dbcon,25);
Q1bw10db = prctile(singlebw10dbcon,25);
Q1bwRMS = prctile(singlebwRMScon,25);
Q1nDur = prctile(singlenDurcon,25);
Q1ppSig = prctile(singleppSignalcon,25);
Q1snr = prctile(singlesnrcon,25);

medpeakFr = prctile(singlepeakFrcon,50);
medcentFr = prctile(singlecentFrcon,50);
meddurClick = prctile(singledurClickcon,50);
medndur95 = prctile(singlendur95con,50);
medndur95Tails = prctile(singlendur95Tailscon,50);
medbw3db = prctile(singlebw3dbcon,50);
medbw10db = prctile(singlebw10dbcon,50);
medbwRMS = prctile(singlebwRMScon,50);
mednDur = prctile(singlenDurcon,50);
medppSig = prctile(singleppSignalcon,50);
medsnr = prctile(singlesnrcon,50);


%Averages 
meandurClick = mean(singledurClickcon);
meannDur = mean(singlenDurcon);
meanndur95 = mean(singlendur95con);
meanndur95Tails = mean(singlendur95Tailscon);
meanbw3db = mean(singlebw3dbcon);
meanbw10db = mean(singlebw10dbcon);
meanbwRMS = mean(singlebwRMScon);
meanpeakFr = mean(singlepeakFrcon);
meancentFr = mean(singlecentFrcon);
meanppSignal = mean(singleppSignalcon);
meansnr = mean(singlesnrcon);

Q3peakFr = prctile(singlepeakFrcon,75);
Q3centFr = prctile(singlecentFrcon,75);
Q3durClick = prctile(singledurClickcon,75);
Q3ndur95 = prctile(singlendur95con,75);
Q3ndur95Tails = prctile(singlendur95Tailscon,75);
Q3bw3db = prctile(singlebw3dbcon,75);
Q3bw10db = prctile(singlebw10dbcon,75);
Q3bwRMS = prctile(singlebwRMScon,75);
Q3nDur = prctile(singlenDurcon,75);
Q3ppSig = prctile(singleppSignalcon,75);
Q3snr = prctile(singlesnrcon,75);


%Standard Deviation
stddurClick = std(singledurClickcon);
stdnDur = std(singlenDurcon);
stdndur95 = std(singlendur95con);
stdndur95Tails = std(singlendur95Tailscon);
stdbw3db = std(singlebw3dbcon);
stdbw10db = std(singlebw10dbcon);
stdbwRMS = std(singlebwRMScon);
stdpeakFr = std(singlepeakFrcon);
stdcentFr = std(singlecentFrcon);
stdppSignal = std(singleppSignalcon);
stdsnr = std(singlesnrcon);

%Make one table, save it as .mat and .xls 
SummaryStats_singleclicks = [Q1durClick Q1ndur95 Q1ndur95Tails Q1bw3db Q1bw10db Q1bwRMS Q1peakFr Q1centFr Q1ppSig Q1snr;...
    meandurClick meanndur95 meanndur95Tails meanbw3db meanbw10db meanbwRMS meanpeakFr meancentFr meanppSignal meansnr;...
    meddurClick medndur95 medndur95Tails medbw3db medbw10db medbwRMS medpeakFr medcentFr medppSig medsnr;...
    Q3durClick Q3ndur95 Q3ndur95Tails Q3bw3db Q3bw10db Q3bwRMS Q3peakFr Q3centFr Q3ppSig Q3snr;...
    stddurClick stdndur95 stdndur95Tails stdbw3db stdbw10db stdbwRMS stdpeakFr stdcentFr stdppSignal stdsnr];


filename = fullfile(inDir,['SummaryStats_singleclicks_',filedate]);
save(filename, 'SummaryStats_singleclicks')
xlswrite([filename,'.xls'],SummaryStats_singleclicks)
%Write csv with headers
headers = {'durClick','ndur95','ndur95Tails','bw3db','bw10db','bwRMS','peakFr','centFr','ppSig','snr'};
csvwrite_with_headers([filename,'.csv'],SummaryStats_singleclicks,headers)

%Make figures for this dataset
%Click duration
figure(1)
hist(singledurClickcon,50)
line([meddurClick meddurClick], [0 25],'Color','r','LineWidth',3);
line([Q1durClick Q1durClick], [0 25],'Color','r','LineWidth',2,'LineStyle',':');
line([Q3durClick Q3durClick], [0 25],'Color','r','LineWidth',2,'LineStyle',':');
text(0.01,0.9,['25th pctile = ',num2str(Q1durClick,3),' \musec'],'Unit','normalized','Color','r')
text(0.01,0.85,['median = ',num2str(meddurClick,3),' \musec'],'Unit','normalized','Color','r')
text(0.01,0.8,['75th pctile = ',num2str(Q3durClick,3),' \musec'],'Unit','normalized','Color','r')
ylim([0 25])
title(['Histogram of single click durations (\musec) (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('Time (\musec)')
filename = fullfile(GraphDir,['SingleClick_Click_duration_',filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
close(figure(1));

%PeakFrequency
figure(3)
vec=0:1:145;
hist(singlepeakFrcon,vec)
%%%Making pretty for publication
xlim([100 145])
set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
set(findobj(gca,'Type','text'),'FontName','Arial','FontSize',12)
set(gca,'FontSize',14,'FontName','Arial')
ylabel('Counts')
xlabel('Frequency (kHz)')
title(['Histogram of single click peak frequencies (kHz)']);

set(gcf,'PaperUnits', 'inches','PaperPosition', [0 0 7.5 10]);
filename = fullfile(GraphDir,['SingleClick_Peak_Frequency_',filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
print('-dpng', filename,'-r600')
close(figure(3));

%Centroid Frequency
figure(3)
vec=0:1:145;
hist(singlecentFrcon,vec)
%%%Making pretty for publication
xlim([100 145])
set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
set(findobj(gca,'Type','text'),'FontName','Arial','FontSize',12)
set(gca,'FontSize',14,'FontName','Arial')
ylabel('Counts')
xlabel('Frequency (kHz)')
title(['Histogram of single click centroid frequencies (kHz)']);

set(gcf,'PaperUnits', 'inches','PaperPosition', [0 0 7.5 10]);
filename = fullfile(GraphDir,['SingleClick_Cent_Frequency_',filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
print('-dpng', filename,'-r600')
close(figure(3));



%ppSignal
figure(4)
hist(singleppSignalcon,50)
line([medppSig medppSig], [0 20],'Color','r','LineWidth',3);
line([Q1ppSig Q1ppSig], [0 20],'Color','r','LineWidth',2,'LineStyle',':');
line([Q3ppSig Q3ppSig], [0 20],'Color','r','LineWidth',2,'LineStyle',':');
text(0.6,0.9,['25th pctile = ',num2str(Q1ppSig,3),' dB'],'Unit','normalized','Color','r')
text(0.6,0.85,['median = ',num2str(medppSig,3),' dB'],'Unit','normalized','Color','r')
text(0.6,0.8,['75th pctile = ',num2str(Q3ppSig,3),' dB'],'Unit','normalized','Color','r')
%ylim([0 300])
%xlim([105 155])
title(['Histgram of single peak-peak amplitudes (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('p-p amplitude (dB)')
filename = fullfile(GraphDir,['SingleClick_P-P_Amplitude_',filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
close(figure(4));

%-3 dB Bandwidth
figure(5)
xbins = 0:1:35;
hist(singlebw3dbcon,xbins)
%xlim([0 20])
text(0.6,0.9,['25th pctile = ',num2str(Q1bw3db,3),' dB'],'Unit','normalized','Color','r')
text(0.6,0.85,['median = ',num2str(medbw3db,3),' dB'],'Unit','normalized','Color','r')
text(0.6,0.8,['75th pctile = ',num2str(Q3bw3db,3),' dB'],'Unit','normalized','Color','r')
title(['Histgram of single -3dB Bandwidths (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('-3dB Bandwidth (kHz)')
% set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
filename = fullfile(GraphDir,['SingleClick_3dBBW_',filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
% saveas(gca, filename, 'fig')
close(figure(5));

%-10 dB Bandwidth
figure(5)
xbins = 0:2:100;
hist(singlebw10dbcon,xbins)
%xlim([0 50])
text(0.6,0.9,['25th pctile = ',num2str(Q1bw10db,3),' dB'],'Unit','normalized','Color','r')
text(0.6,0.85,['median = ',num2str(medbw10db,3),' dB'],'Unit','normalized','Color','r')
text(0.6,0.8,['75th pctile = ',num2str(Q3bw10db,3),' dB'],'Unit','normalized','Color','r')
title(['Histgram of single -10dB Bandwidths (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('-10dB Bandwidth (kHz)')
% set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
filename = fullfile(GraphDir,['SingleClick_10dBBW_',filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
% saveas(gca, filename, 'fig')
close(figure(5));


%RMS Bandwidth
figure(6)
xbins = 0:1:30;
hist(singlebwRMScon,xbins)
xlim([0 30])
% ylim([0 100])
text(0.6,0.9,['25th pctile = ',num2str(Q1bwRMS,3),' dB'],'Unit','normalized','Color','r')
text(0.6,0.85,['median = ',num2str(medbwRMS,3),' dB'],'Unit','normalized','Color','r')
text(0.6,0.8,['75th pctile = ',num2str(Q3bwRMS,3),' dB'],'Unit','normalized','Color','r')
title(['Histgram of single RMS Bandwidths (n = ',strnclicks,')']);
ylabel('Counts')
xlabel('RMS Bandwidth (kHz)')
% set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
filename = fullfile(GraphDir,['SingleClick_RMSBW_',filedate]);
saveas(gca, filename, 'tif')
saveas(gca, filename, 'jpg')
% saveas(gca, filename, 'fig')
close(figure(6));






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Then combine single spectra and calculate one mean, plot with the mean noise.
%First, calculate the mean. start by looping through all the entries in the
%cell array and saving those to a non-cell matrix
numspecClick = size(singlespecClickTfcon,1);
numspecNoise = size(singlespecNoiseTfcon,1);
concatspecs = [];
concatspecsN = [];
for s = 1:numspecClick
        concatspecs = [concatspecs;singlespecClickTfcon{s,1}'];
end
grandmeanSpec = mean(concatspecs);

for s = 1:numspecNoise
        concatspecsN = [concatspecsN;singlespecNoiseTfcon{s,1}'];
end
%Get unique noise spectra (because some/many are repeats)
concatspecsNU = unique(concatspecsN,'rows');
grandmeanSpecN = mean(concatspecsNU);

%Make frequency axis
freqs = f;

% figure(8)
% plot(freqs,grandmeanSpec(1,1:end), 'Color','k','LineWidth',1);
% hold on
% plot(freqs,grandmeanSpecN(1,1:end), 'Color',[0.7, 0.7, 0.7], 'LineWidth',1);
% % plot(freqs,singlemeanSpecNoises(1,2:end), 'Color',[0.7, 0.7, 0.7], 'LineWidth',1);
% xlim([0 175]);
% legend('mean click','mean noise','Location','southwest')
% title(['Mean Spectrum of single Clicks (n = ',num2str(numspecClick),')']);
% xlabel('Frequency (kHz)')
% ylabel('Amplitude (dB)')
% filename = fullfile(GraphDir,['SingleClick_MeanSpecGrand',filedate]);
% saveas(gca, filename, 'tif')
% saveas(gca, filename, 'jpg')
% close(figure(8));

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
toohighfreq = find(plotfreqs > 170); 
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

filename = fullfile(GraphDir,['SingleClick_MeanSpecGrand_forPub',filedate]);
print('-dtiff',filename)
%print('-djpeg', filename)
print('-dpng', filename,'-r600')
%print('-deps', filename,'-r200')
% saveas(gca, filename, 'fig')
close(figure(9));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

