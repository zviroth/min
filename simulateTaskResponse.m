close all
clear all
upsampleFactor = 10;
TR = 1.5;
%set canonical HRF
modelParams = struct;
sampleDuration = TR/upsampleFactor;
sampleDelay=sampleDuration/2;
defaultParams=1;
modelParams.z = 10;
[modelParams hrfModel] = hrfDoubleGamma(modelParams,sampleDuration,sampleDelay,defaultParams);
i=0;
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
numTrials = 15;
trialLength = 10;
runsPerRwd = 100;

% snr = [1 0.95];
% temporalJitter = [0 1.5];%std of jitter, in Fourier phase

snr = [0.1:0.1:3];
temporalJitter = zeros(size(snr));

% temporalJitter = [0 : 0.01: 0.1];%std of jitter, in Fourier phase
% snr = 1*ones(size(temporalJitter));


numRwd=length(snr);
for rwd=1:numRwd
    plotColors{rwd} = [1 - (rwd-1)/(numRwd-1), 0, (rwd-1)/(numRwd-1)];
end
% snr = [1 1];
% temporalJitter = [0 1.5];%std of jitter, in Fourier phase



taskAmp = 1;%ones(numRwd,1);
T = numTrials*trialLength;
upT = T*upsampleFactor;
for rwd=1:numRwd
    for r=1:runsPerRwd
        runTC = zeros(1,upT);
        taskTiming = 1:trialLength*upsampleFactor:upT;%temporal jitter will be added here, after upsampling
        noisyTiming = taskTiming + (temporalJitter(rwd)*upsampleFactor*2*pi/trialLength)*randn(size(taskTiming));
        noisyTiming(1) = 1;
        
        runTC(ceil(noisyTiming)) = ones;
        temp = conv(runTC,hrfModel);
        runTC = temp(1:upT);%crop end
        
        rwdSignal(rwd,:,r) = taskAmp*runTC(1:upsampleFactor:end);%downsample
%         rwdSignal(rwd,:,r) = filterTC(runTC(1:upsampleFactor:end));%downsample

    end
    n = (taskAmp/snr(rwd))*randn(size(rwdSignal(rwd,:,:)));
    oneOverF = (log10(1:T))*1./(1:T);
    rwdNoise(rwd,:,:) = abs(ifft(repmat(oneOverF,1,1,runsPerRwd).*fft(n)));
    
end

%%
rwdTC = rwdSignal + rwdNoise;
rwdTC = zscore(rwdTC,0,2);

rwdTrials = reshape(rwdTC,numRwd,trialLength,[]);
meanTrial = mean(rwdTrials,3);
meanRun = squeeze(mean(rwdTC,3));
for rwd=1:numRwd
    f = squeeze(fft(rwdTrials(rwd,:,:)));
    fftTrialAmp(rwd,:) = abs(f(2,:));
    fftTrialPh(rwd,:) = angle(f(2,:));
    fftRun(rwd,:) = abs(fft(meanRun(rwd,:)));
end

nframes = length(fftRun(rwd,:));
%from plotMeanFourierAmp.m
fftRun = fftRun / (nframes/2);
frequencies = [0:nframes-1]/(nframes*TR);
fftRun = fftRun(:,1:floor(nframes/2));
frequencies = frequencies(1:floor(nframes/2));


trialsAmp = squeeze(std(rwdTrials,0,2));
ampMeanTrial = std(meanTrial,0,2);
meanTrialsAmp = mean(trialsAmp,2);

phVar = circ_std(fftTrialPh');
ampVar = std(fftTrialAmp,0,2);
totalVar = mean(std(rwdTrials,0,3),2);



%% plot STD of mean trial and mean of trial STD as function of SNR
i=i+1; figure(i) ;clf
rows=2;
cols=4;
subplot(rows,cols,1)
plot(snr, ampMeanTrial)
legend('amplitude of mean trial');
subplot(rows,cols,2)
plot(snr, meanTrialsAmp)
legend('mean trials amplitude');
subplot(rows,cols,3)
plot(snr, phVar)
legend('FFT phase variability');
subplot(rows,cols,4)
plot(snr, totalVar)
legend('mean timepoint variability');

subplot(rows,cols,cols+1)
plot(temporalJitter, ampMeanTrial)
legend('amplitude of mean trial');
subplot(rows,cols,cols+2)
plot(temporalJitter, meanTrialsAmp)
legend('mean trials amplitude');
subplot(rows,cols,cols+3)
plot(temporalJitter, phVar)
legend('FFT phase variability');
subplot(rows,cols,cols+4)
plot(temporalJitter, totalVar)
legend('mean timepoint variability');


%% plot data
i=i+1; figure(i) ;clf
rows=2;
cols=numRwd;
for rwd=1:numRwd
    subplot(rows,cols,rwd)
    plot(squeeze(mean(rwdTC(rwd,:,:),3)))
    subplot(rows,cols,rwd+cols)
    plot(meanTrial(rwd,:))
end
ampMeanTrial;
meanTrialsAmp;

%% 
linewidth = 1;
i=i+1; figure(i) ;clf
rows=1; cols=2;
subplot(rows,cols,1)
for rwd=1:numRwd
    plot(meanTrial(rwd,:),'color',plotColors{rwd},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,2)
for rwd=1:numRwd
    plot(frequencies,fftRun(rwd,:),'.-','color',plotColors{rwd},'linewidth',linewidth);
    hold on
end


%%
function filteredTC = filterTC(tc)
filt = ones(1,length(tc));
filtNonzeros = [0,0,0.117503097415405,0.588887709492812,0.904365555167461,0.988891003461758,0.999355620179431];
filt(2:8) = filtNonzeros;
filt(end-6:end) = filtNonzeros;
filteredTC = real(ifft(fft(tc) .* filt ))';
end