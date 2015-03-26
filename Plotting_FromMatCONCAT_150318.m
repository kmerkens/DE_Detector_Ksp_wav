%Plotting when you just load a single .mat file into the workspace. 

% fs = 480000;
% fileID = 'dalls_20080802_080000_7_30_00';
% GraphDir = 'C:\Users\Karlina.Merkens\Documents\Porpoise\OtherRecordings\TYack_Dalls_Wild\metadata\dalls';

% fs = 500000;
% fileID = 'kogiasima_20020716_170033_lab';
% GraphDir = 'C:\Users\Karlina.Merkens\Documents\Kogia\OtherRecordings\DMann_Ksima_captive\metadata\kogia';

fs = 375000;
fileID = 'ClicksOnlyConcat';
GraphDir = 'C:\Users\Karlina.Merkens\Documents\Kogia\OtherRecordings\VJanik_Ksima_Wild\metadata\kogia';



%Calculate ici
pos1 = [clickDnum(:,1);0];
pos2 = [0;clickDnum(:,1)];
ici = pos1(2:end-1)-pos2(2:end-1);
iciVec = datevec(ici); %To get the vector
iciSecs = iciVec(:,6);
       
 %If it's bigger than 0.5 seconds, remove it - there's a gap between click
 %bunches
iciTooBig = find(iciSecs > 0.25);
iciSecs(iciTooBig) = [];
icims = iciSecs*1000;

params = {'median peak frequency (kHz)','median ipi (ms)',...
    'median duration (us)','median received level (dB re 1 uPa)',...
    'median center frequency (kHz)'};
percentiles = [10 50 90];
medianValue(1) = prctile(peakFrcon,50);%calculate median peak frequency
medianValue(2) = prctile(icims,50);%calculate median inter-pulse interval
medianValue(3) = prctile(durClickcon,50);%calculate median duration
medianValue(4) = prctile(ppSignalcon,50);%calculate median inter-pulse interval
%medianValue(5) = prctile(F0Sel,50);%calculate median center frequency
%medianValue(5) = encStart;
clickCount = sum(clickDnum);%count number of clicks in analysis
maxRL = max(ppSignalcon);

%sort spectras for peak frequency and prepare for plotting concetanated
%spectrogram
[a b]=sort(peakFrcon);

specSorted=[];
for c=1:length(b)
    thisspec = cell2mat(specClickTfcon(b(c),:));
    specSorted(c,:)=thisspec;
end
specSorted=specSorted.';

N=size(specSorted,1)*2;
%f=0:(fs/2000)/(N/2-1):fs/2000;
datarow=size(specSorted,2);

%calculate mean spectra for click and noise
meanSpecClick = mean(specSorted');

%meanSpecNoise=mean(specNoiseSel);

%         sep = strfind(pathstr,'\');
%         disk = pathstr(sep(2)+1:length(pathstr));
% 
%         figure('Name', sprintf('%s %s', disk, datestr(encStart)),...
%              'Position',([0,0,1200,800]))

subplot(2, 2, 1);
vec=0:1:160;
hist(peakFrcon,vec)
xlim([0 200])
%ylim([0 4])
xlabel('peak frequency (kHz)')
ylabel('counts')
text(0.05,0.9,['pfr =',num2str(medianValue(1)),' kHz'],'Unit','normalized')
%text(0.5,0.8,['cfr =',num2str(medianValue(5)),' kHz'],'Unit','normalized')

subplot(2,2,2)
vec=0:5:260;
hist(icims,vec)
xlim([0 260])
%ylim([0 7])
xlabel('inter-pulse interval (ms)')
ylabel('counts')
text(0.5,0.9,['dur =',num2str(medianValue(3)),' \mus'],'Unit','normalized')
text(0.5,0.8,['ipi =',num2str(medianValue(2)),' ms'],'Unit','normalized')

subplot(2,2,3)
plot(f,meanSpecClick,'LineWidth',2), hold on
%plot(f,meanSpecNoise,':k','LineWidth',2), hold off
xlabel('Frequency (kHz)'), ylabel('Normalized amplitude (dB)')
ylim([0 60])
xlim([0 200])
title(['Mean click spectra, n=',num2str(size(specSorted,2))],'FontWeight','bold')
text(0.5,0.9,['ppRL =',num2str(medianValue(4))],'Unit','normalized')

subplot(2,2,4)
imagesc(1:datarow, f, specSorted); axis xy; colormap(gray)
xlabel('Click number'), ylabel('Frequency (kHz)')
title(['Clicks sorted by peak frequency'],'FontWeight','bold')

seqIdent = ([fileID,'_multi']);
filename = fullfile(GraphDir,seqIdent);
saveas(gca, filename, 'jpg')

close

figure
plot(f,meanSpecClick,'LineWidth',2), hold on
%plot(f,meanSpecNoise,':k','LineWidth',2), hold off
xlabel('Frequency (kHz)'), ylabel('Normalized amplitude (dB)')
ylim([0 60])
xlim([0 200])
title(['Mean click spectra, n=',num2str(size(specSorted,2))],'FontWeight','bold')
text(0.05,0.9,['pfr =',num2str(medianValue(1)),' kHz'],'Unit','normalized')

seqIdent = ([fileID,'_spec']);
filename = fullfile(GraphDir,seqIdent);
saveas(gca, filename, 'jpg')
saveas(gca, filename, 'tif')

close
