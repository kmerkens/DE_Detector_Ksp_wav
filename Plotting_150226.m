%Plotting

figure(1)
hist(peakFrcon, 20)
title('peak Frequency')

figure(2)
hist(durClickcon,20)
title('Click duration')

figure(3)
hist(ppSignalcon,20)
title('peak to peak signal')

%Make a concat spectrogram of all the specClickTfcon cells