close all
clear all
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
end
zScoreString = '';
if toZscore
    zScoreString = '_zscored';
end
globalMeanString = '';
if regressGlobalMean
    globalMeanString = '_globalRegressed';
end
ConcatProjStr = '';
if ConcatProj
    ConcatProjStr = 'ConcatProj';
end
load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], ...
    'subFolders', 'subTrialResponse', 'subResponse');

upsampleFactor = 10;

goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% goodSubs = [1:3];
iRoi=2;
for iSub = 1:length(goodSubs)%length(subdirs)
    for rwd=1:2
        reshapedTrials = subTrialResponse{goodSubs(iSub),iRoi,rwd};
        %trial-by-trial variability per subject
        subTimepointVar(iSub,rwd,:) = std(reshapedTrials,0,2);
        subRwdHrf(iSub,rwd,:) = subResponse(goodSubs(iSub),iRoi,rwd,:);
    end
    subHrf(iSub,:) = 0.5*(subRwdHrf(iSub,1,:)+subRwdHrf(iSub,2,:));
    temp = interp1(1:length(subHrf(iSub,:))+1,[subHrf(iSub,:) subHrf(iSub,1)],1:(1/upsampleFactor):(length(subHrf(iSub,:))+1));
    subInterpHrf(iSub,:) = temp(1:end-1);
    subInterpHrf(iSub,:) = subInterpHrf(iSub,:) - subInterpHrf(iSub,1);%make it start at 0

end


toZscore=0;
TR = 1.5;
%set canonical HRF
modelParams = struct;
sampleDuration = TR/upsampleFactor;
sampleDelay=sampleDuration/2;
defaultParams=1;
% params that make the variability timecourse similar to 
modelParams.x = 6;%9;
modelParams.y = 16;
modelParams.z = 6;
% params that make the variability timecourse similar to temporal jitter
modelParams.x = 12;
modelParams.y = 16;
modelParams.z = 6;

% % no negative dip
% modelParams.x = 6;%9;
% modelParams.y = 16;
% modelParams.z = 100;

[modelParams hrfModel] = hrfDoubleGamma(modelParams,sampleDuration,sampleDelay,defaultParams);
i=0;
% plot(hrfModel);
% keyboard

%%
plotColorsSnr = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
numTrials = 16;
trialLength = 10;
runsPerRwd = 10;
T = numTrials*trialLength;
upT = T*upsampleFactor;
oneOverF(1,1,1,1,:) = (log10(1:T))*1./(1:T);
oneOverF(1,1,1,1,2:end) = oneOverF(1,1,1,1,1:end-1);%we only care about the second component onwards
oneOverF(1,1,1,1,T/2+1:end) = oneOverF(1,1,1,1,1+T/2:-1:2);%make it symmetric


%%

npoints = 20;
snr = logspace(2, -0.5, npoints);
snr(1) = Inf;
snr(3:end) = [];
jitter = linspace(0, pi/2, npoints);
% jitter = jitter(1:9);
ampJitter = logspace(-2,0,npoints);
ampJitter(1) = 0;


for isnr=1:length(snr)
    plotColorsSnr{isnr} = [1 - (isnr-1)/(length(snr)-1), 0, (isnr-1)/(length(snr)-1)];
end
for ijitter=1:length(jitter)
    plotColorsJitter{ijitter} = [1 - (ijitter-1)/(length(jitter)-1), 0, (ijitter-1)/(length(jitter)-1)];
end
for iampJitter=1:length(ampJitter)
    plotColorsAmpJitter{iampJitter} = [1 - (iampJitter-1)/(length(ampJitter)-1), 0, (iampJitter-1)/(length(ampJitter)-1)];
end
% snr = [1 1];
% temporalJitter = [0 1.5];%std of jitter, in Fourier phase



taskAmp = 1;%ones(numRwd,1);


for isnr=1:length(snr)
    for ijitter = 1:length(jitter)
        for iampJitter=1:length(ampJitter)
            for r=1:runsPerRwd
                runTC = zeros(1,upT);
                taskTiming = 1:trialLength*upsampleFactor:upT;%temporal jitter will be added here, after upsampling
                noisyTiming = taskTiming + (jitter(ijitter)*upsampleFactor*trialLength/(2*pi))*randn(size(taskTiming));
                noisyTiming(noisyTiming<1) = 1;
                
                noisyTiming(noisyTiming>T*upsampleFactor) = T;
                
                runTC(ceil(noisyTiming)) = ones;
%                 runTC(ceil(noisyTiming)) = taskAmp*(1+atan(ampJitter(iampJitter)*randn(numTrials,1))./(pi/2));
                runTC(ceil(noisyTiming)) = taskAmp*(1+atan(ampJitter(iampJitter)*randn(numTrials,1)));
                
                for iSub=1:length(goodSubs)
%                     hrfModel
                    temp = conv(runTC,subInterpHrf(iSub,:));
%                     temp = conv(runTC,hrfModel);
                    convRunTC = temp(1:upT);%crop end
                    rwdSignal(iSub,isnr,ijitter,iampJitter,:,r) = convRunTC(1:upsampleFactor:end);%downsample 
                end

            end
            n = (taskAmp/snr(isnr))*randn(size(rwdSignal(iSub,isnr,ijitter,iampJitter,:,:)));
            
            %WITHOUT ABS
            temp = ifft(repmat(oneOverF,1,1,1,1,1,runsPerRwd).*fft(n));
            rwdNoise(:,isnr,ijitter,iampJitter,:,:) = repmat(temp,length(goodSubs),1,1,1,1,1);%same for all subjects
        end
    end
end

%%
rwdTC = rwdSignal + rwdNoise;
rwdTC = rwdTC*upsampleFactor;
if toZscore
    rwdTC = zscore(rwdTC,0,5);
end
rwdTC = rwdTC(:,:,:,:,trialLength+1:end,:);%junk first cycle
for iSub=1:length(goodSubs)
    for isnr=1:length(snr)
        for ijitter = 1:length(jitter)
            for iampJitter=1:length(ampJitter)
                for r=1:runsPerRwd
                    rwdTC(iSub,isnr,ijitter,iampJitter,:,r) = filterTC(squeeze(rwdTC(iSub,isnr,ijitter,iampJitter,:,r)));
                end
            end
        end
    end
end


% for iSub=1:length(goodSubs)
%     rwdTrials(iSub,:,:,:,:,:) = reshape(rwdTC(iSub,:,:,:,:,:),length(snr),length(jitter),length(ampJitter),trialLength,[]);
% end
rwdTrials = reshape(rwdTC,length(goodSubs),length(snr),length(jitter),length(ampJitter),trialLength,[]);%(isnr,ijitter,iampjitter,trialLength,numTrials)
meanTrial = mean(rwdTrials,6);
meanRun = squeeze(mean(rwdTC,6));
for iSub=1:length(goodSubs)
    for isnr=1:length(snr)
        for ijitter = 1:length(jitter)
            for iampJitter=1:length(ampJitter)
                f = squeeze(fft(rwdTrials(iSub,isnr,ijitter,iampJitter,:,:)));
                fftTrialAmp(iSub,isnr,ijitter,iampJitter,:) = abs(f(2,:));
                fftTrialPh(iSub,isnr,ijitter,iampJitter,:) = angle(f(2,:));
                fftRun(iSub,isnr,ijitter,iampJitter,:) = abs(fft(meanRun(iSub,isnr,ijitter,iampJitter,:)));
            end
        end
    end
end

nframes = length(fftRun);
%from plotMeanFourierAmp.m
fftRun = fftRun / (nframes/2);
frequencies = [0:nframes-1]/(nframes*TR);
fftRun = fftRun(:,:,:,:,1:floor(nframes/2));
frequencies = frequencies(1:floor(nframes/2));



%%
trialsAmp = squeeze(std(rwdTrials,0,5));
ampMeanTrial = std(meanTrial,0,5);

fftMeanTrial = fft(meanTrial,[],5);
f = fftMeanTrial(:,:,:,:,2);
phMeanTrial = angle(f);


meanTrialsAmp = mean(trialsAmp,5);
stdTrialsAmp = std(trialsAmp,0,5);

meanTrialsFFTamp = mean(fftTrialAmp,5);

phVar = circ_std(fftTrialPh,[],[],5);
ampVar = std(fftTrialAmp,0,5);
totalVar = mean(std(rwdTrials,0,6),5);

%% timepoint variability timecourse
% [simMeanAmpTimepointVar simMeanTempTimepointVar];

simTimepointVar = squeeze(std(rwdTrials,0,6));
simMeanVar = mean(simTimepointVar(:,1,:,:,:),5);
simMeanAmpTimepointVar = squeeze(mean(simMeanVar(:,1,1,:)));
simMeanTempTimepointVar = squeeze(mean(simMeanVar(:,1,:,1)));
% ampNoiseLevels = [8 10; 9 10];
% tempNoiseLevels = [3 5; 4 5];
ampNoiseLevels = [ 19 20; 17 20;  13 17; 13 20];
tempNoiseLevels = [ 7 9; 5 9; 3 5; 3 9];
% rwdNoiseLevel = [1 10; 2 9; 3 8; 4 7; 5 6];
% rwdNoiseLevel = [1 2; 3 4; 5 6; 7 8;  9 10];
% rwdNoiseLevel = [2 3; 4 5; 6 7; 8 9];

rows=size(ampNoiseLevels,1);
cols=5;
i=0;
clear pval
for inoiseContrast = 1:size(ampNoiseLevels,1)
    %amplitude variability
    isnr=1; ijitter=1;
    simAmpTimepointVar = squeeze(std(rwdTrials(:,isnr,ijitter,ampNoiseLevels(inoiseContrast,:),:,:),0,6));
    simAmpTimepointVarDiff = squeeze(simAmpTimepointVar(:,2,:) - simAmpTimepointVar(:,1,:));
    subTimepointVarDiff = squeeze(subTimepointVar(:,2,:) - subTimepointVar(:,1,:));
    ampVarCorr = diag(corr(simAmpTimepointVarDiff', subTimepointVarDiff'));
    
    %temporal variability
    isnr=1; iampJitter=1;
    simTempTimepointVar = squeeze(std(rwdTrials(:,isnr,tempNoiseLevels(inoiseContrast,:),iampJitter,:,:),0,6));
    simTempTimepointVarDiff = squeeze(simTempTimepointVar(:,2,:) - simTempTimepointVar(:,1,:));
    tempVarCorr = diag(corr(simTempTimepointVarDiff', subTimepointVarDiff'));
    
    [h pval(inoiseContrast)] = ttest(atan(tempVarCorr) - atan(ampVarCorr));
    
    i=i+1;
    subplot(rows,cols,i)
    % plot(subHrf')
    plot(subInterpHrf')
    title('subject HRF');
    i=i+1;
    subplot(rows,cols,i)
    plot(subTimepointVarDiff')
    title('subject variability');
    i=i+1;
    subplot(rows,cols,i)
    plot(simAmpTimepointVarDiff')
    plot(squeeze(simAmpTimepointVar(:,2,:))')
    title('timepoint variability - amplitude jitter');
    i=i+1;
    subplot(rows,cols,i)
    plot(simTempTimepointVarDiff')
    plot(squeeze(simTempTimepointVar(:,2,:))')
    title('timepoint variability - temporal jitter');
    i=i+1;
    subplot(rows,cols,i)
    plot(-1:1:1,-1:1:1,'-'); hold on
    plot(ampVarCorr,tempVarCorr,'.','markersize',25);
    xlabel('amplitude jitter');
    ylabel('temporal jitter');
%     title([num2str(rwdNoiseLevel(inoiseContrast,2)) ' - ' num2str(rwdNoiseLevel(inoiseContrast,1))]);
    title([num2str(ampNoiseLevels(inoiseContrast,2)) ' - ' num2str(ampNoiseLevels(inoiseContrast,1)) ', '...
        num2str(tempNoiseLevels(inoiseContrast,2)) ' - ' num2str(tempNoiseLevels(inoiseContrast,1))]);
    
end
pval
% set(gcf,'position',[200 200 1400 300*rows]);
set(gcf,'position',[200 200 1400 1400]);

%%
function filteredTC = filterTC(tc)
filt = ones(1,length(tc));
filtNonzeros = [0,0,0.117503097415405,0.588887709492812,0.904365555167461,0.988891003461758,0.999355620179431];
filt(2:8) = filtNonzeros;
filt(end-6:end) = filtNonzeros;
filteredTC = real(ifft(fft(tc) .* filt' ))';

% filteredTC = tc;
end