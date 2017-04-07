noisetest = bpDataHi(1:313);
clicktest = bpDataHi(496:555);

nsamp = length(clicktest);

noisetest = bpDataHi(1:nsamp); %Recalculate noisetest

noisetestcts = noisetest*2^15;
clicktestcts = clicktest*2^15;

clickwinLength = length(clicktestcts);
noisewinLength = length(noisetestcts);

clickwind = hann(clickwinLength);
noisewind = hann(noisewinLength);

wClick = clicktestcts.*clickwind';
wNoise = noisetestcts.*noisewind';  

fftSize = ceil(hdr.fs * p.frameLengthUs / 1E6);

spClick = 20*log10(abs(fft(wClick,fftSize)));
spNoise = 20*log10(abs(fft(wNoise,fftSize)));

sub = 10*log10(hdr.fs/fftSize);

spClickSub = spClick-sub;
spNoiseSub = spNoise-sub; 

spClickSub = spClickSub(:,1:fftSize/2);
spNoiseSub = spNoiseSub(:,1:fftSize/2);

plot(spClickSub)
hold on
plot(spNoiseSub,'r')
