close all
clear all
upsampleFactor = 10;
TR = 1.5;
%set canonical HRF
modelParams = struct;
sampleDuration = TR/upsampleFactor;
sampleDelay=sampleDuration/2;
defaultParams=1;
% modelParams.x = 10;
% modelParams.z = 6;
[modelParams hrfModel] = hrfDoubleGamma(modelParams,sampleDuration,sampleDelay,defaultParams);
i=0;


plotColorsSnr = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
numTrials = 16;
trialLength = 10;
runsPerRwd = 100;
T = numTrials*trialLength;
upT = T*upsampleFactor;
oneOverF(1,1,:) = (log10(1:T))*1./(1:T);
oneOverF(1,1,2:end) = oneOverF(1,1,1:end-1);%we only care about the second component onwards
oneOverF(1,1,T/2+1:end) = oneOverF(1,1,1+T/2:-1:2);%make it symmetric


i=i+1; figure(i) ;clf
subplot(1,2,1)
plot(hrfModel(1:upsampleFactor:end))
hold on
plot(abs(diff(hrfModel(1:upsampleFactor:end))))
vline(10);

subplot(1,2,2)
plot(squeeze(oneOverF))


%%

% snr = [1 0.95];
% temporalJitter = [0 1.5];%std of jitter, in Fourier phase
npoints = 10;
% nsr = linspace(0.01, 5, npoints);
% snr = 1./nsr;
snr = linspace(5, 0.01, npoints);
% snr = [3:-0.1:0.1];
% jitter = zeros(size(snr));

jitter = linspace(0, pi/2, npoints);
% jitter = [0 : 0.2: 2*pi];%std of jitter, in Fourier phase
% snr = 1*ones(size(temporalJitter));



for isnr=1:length(snr)
    plotColorsSnr{isnr} = [1 - (isnr-1)/(length(snr)-1), 0, (isnr-1)/(length(snr)-1)];
end

for ijitter=1:length(jitter)
    plotColorsJitter{ijitter} = [1 - (ijitter-1)/(length(jitter)-1), 0, (ijitter-1)/(length(jitter)-1)];
end
% snr = [1 1];
% temporalJitter = [0 1.5];%std of jitter, in Fourier phase



taskAmp = 1;%ones(numRwd,1);


for isnr=1:length(snr)
    for ijitter = 1:length(jitter)
        for r=1:runsPerRwd
            runTC = zeros(1,upT);
            taskTiming = 1:trialLength*upsampleFactor:upT;%temporal jitter will be added here, after upsampling
            noisyTiming = taskTiming + (jitter(ijitter)*upsampleFactor*trialLength/(2*pi))*randn(size(taskTiming));
            noisyTiming(noisyTiming<1) = 1;
            noisyTiming(noisyTiming>T*upsampleFactor) = T;
            
            runTC(ceil(noisyTiming)) = ones;
            temp = conv(runTC,hrfModel);
            runTC = temp(1:upT);%crop end
            
            rwdSignal(isnr,ijitter,:,r) = taskAmp*runTC(1:upsampleFactor:end);%downsample
            
        end
        n = (taskAmp/snr(isnr))*randn(size(rwdSignal(isnr,ijitter,:,:)));
        
        rwdNoise(isnr,ijitter,:,:) = abs(ifft(repmat(oneOverF,1,1,1,runsPerRwd).*fft(n)));
    end
end

%%
rwdTC = rwdSignal + rwdNoise;
rwdTC = zscore(rwdTC,0,3);
rwdTC = rwdTC(:,:,trialLength+1:end,:);%junk first cycle
for isnr=1:length(snr)
    for ijitter = 1:length(jitter)
        for r=1:runsPerRwd
            rwdTC(isnr,ijitter,:,r) = filterTC(squeeze(rwdTC(isnr,ijitter,:,r)));
        end
    end
end



rwdTrials = reshape(rwdTC,length(snr),length(jitter),trialLength,[]);
meanTrial = mean(rwdTrials,4);
meanRun = squeeze(mean(rwdTC,4));
for isnr=1:length(snr)
    for ijitter = 1:length(jitter)
        f = squeeze(fft(rwdTrials(isnr,ijitter,:,:)));
        fftTrialAmp(isnr,ijitter,:) = abs(f(2,:));
        fftTrialPh(isnr,ijitter,:) = angle(f(2,:));
        fftRun(isnr,ijitter,:) = abs(fft(meanRun(isnr,ijitter,:)));
    end
end

nframes = length(fftRun);
%from plotMeanFourierAmp.m
fftRun = fftRun / (nframes/2);
frequencies = [0:nframes-1]/(nframes*TR);
fftRun = fftRun(:,:,1:floor(nframes/2));
frequencies = frequencies(1:floor(nframes/2));


trialsAmp = squeeze(std(rwdTrials,0,3));
ampMeanTrial = std(meanTrial,0,3);
meanTrialsAmp = mean(trialsAmp,3);


phVar = circ_std(fftTrialPh,[],[],3);
ampVar = std(fftTrialAmp,0,3);
totalVar = mean(std(rwdTrials,0,4),3);



%% plot STD of mean trial and mean of trial STD as function of SNR
i=i+1; figure(i) ;clf
rows=1;
cols=4;
subplot(rows,cols,1)
imagesc(ampMeanTrial)
title('amplitude of mean trial');
subplot(rows,cols,2)
imagesc(meanTrialsAmp)
title('mean trials amplitude');
subplot(rows,cols,3)
imagesc(phVar)
title('FFT phase variability');
subplot(rows,cols,4)
imagesc(squeeze(totalVar))
title('mean timepoint variability');

for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    xlabel('temporal jitter (radians)');
    xticks(1:10:length(jitter))
    xticklabels(jitter(1:10:end));
    axis('image')
    yticks([])
end
subplot(rows,cols,1)
ylabel('SNR');
yticks(1:5:length(snr))
yticklabels(snr(1:5:end));
% %% plot data
% i=i+1; figure(i) ;clf
% rows=2;
% cols=numRwd;
% for rwd=1:numRwd
%     subplot(rows,cols,rwd)
%     plot(squeeze(mean(rwdTC(rwd,:,:),3)))
%     subplot(rows,cols,rwd+cols)
%     plot(meanTrial(rwd,:))
% end
% ampMeanTrial;
% meanTrialsAmp;
% 

%% 
linewidth = 0.1;
markersize = 10;
i=i+1; figure(i) ;clf
rows=2; cols=3;
subplot(rows,cols,1)
ijitter=1;
for isnr=1:length(snr)    
    plot(squeeze(meanTrial(isnr,ijitter,:)),'color',plotColorsSnr{isnr},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,2)
for isnr=1:length(snr)
    plot(frequencies,squeeze(fftRun(isnr,ijitter,:)),'.-','color',plotColorsSnr{isnr},'linewidth',linewidth,'markersize',markersize);
%     plot(frequencies,log(squeeze(fftRun(isnr,ijitter,:))),'.-','color',plotColorsSnr{isnr},'linewidth',linewidth,'markersize',markersize);
    hold on
end
ylabel('log(FFT amplitude)');
subplot(rows,cols,3)
for isnr=1:length(snr)
    plot(squeeze(std(rwdTrials(isnr,ijitter,:,:),0,4)),'color',plotColorsSnr{isnr},'linewidth',linewidth);
    hold on
end

subplot(rows,cols,cols+1)
isnr=1;
for ijitter=1:length(jitter)    
    plot(squeeze(meanTrial(isnr,ijitter,:)),'color',plotColorsJitter{ijitter},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,cols+2)
for ijitter=1:length(jitter)
    plot(frequencies,squeeze(fftRun(isnr,ijitter,:)),'.-','color',plotColorsJitter{ijitter},'linewidth',linewidth,'markersize',markersize);
%     plot(frequencies,log(squeeze(fftRun(isnr,ijitter,:))),'.-','color',plotColorsJitter{ijitter},'linewidth',linewidth,'markersize',markersize);
    hold on
end
ylabel('log(FFT amplitude)');
subplot(rows,cols,cols+3)
for ijitter=1:length(jitter)
    plot(squeeze(std(rwdTrials(isnr,ijitter,:,:),0,4)),'color',plotColorsSnr{ijitter},'linewidth',linewidth);
    hold on
end
%%
i=i+1; figure(i) ;clf
rows=1;
cols=2;
subplot(rows,cols,1)
plot(squeeze(rwdTC(1,1,:,1)))
subplot(rows,cols,2)
plot(squeeze(rwdTC(1,end,:,1)))

%%
function filteredTC = filterTC(tc)
filt = ones(1,length(tc));
filtNonzeros = [0,0,0.117503097415405,0.588887709492812,0.904365555167461,0.988891003461758,0.999355620179431];
filt(2:8) = filtNonzeros;
filt(end-6:end) = filtNonzeros;
filteredTC = real(ifft(fft(tc) .* filt' ))';

end