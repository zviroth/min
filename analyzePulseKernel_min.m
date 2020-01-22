close all
mrQuit
clear all
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';

load([dataFolder 'rwd_physioRegress.mat'], 'dataFolder', 'subFolders', 'numSubs','roiNames', ...
    'numRuns','numTRs','concatInfo',...
    'trialsPerRun','trialLength','junkedFrames','TR',...
    'ecgselect','ecgSampleRate','ecgTrial','ecgRunLength','ecgInterpMethod',...
    'numRuns','concatInfo','numTRs','globalMean',...
    'ecg','ecgPeaksDiff','ecgPeaksAmp','ecgRateTime','interpPeaksDiff',...
    'trialsPeakDiff','ecgPulseRate','interpPulseRate',...
    'trialsPulseRate','meanTrialPulseRateRun','meanTrialPeaksDiffRun',...
    'rwdMeanPeaksDiff','rwdMeanPulseRate','rwdPulseTC','downsampledPulse',...
    'designMatPulse',...
    'origVar','origVarRegressGlobal','pulseResidualVar','pulseResidualVarRegressGlobal',...
    'meanRoiTC','meanRoiTCregress',...
    'pulseRsqr','pulseRsqrRegressGlobal',...
    'pulseKernel','pulseKernelRegress',...
    'pulseResidualTC', 'pulseResidualTCregressGlobal',...
    'rv','rvt','respselect','resp','respPeaks','respPeaksDiff',...
    'respPeaksAmp','interpPeaksAmp','respRateTime',...
    'interpRespPeaksDiff','respTroughs','respTroughsAmp','interpTroughsAmp',...
    'physioKernel','physioKernelRegress','physioResidualTC','physioResidualVar',...
    'physioResidualTCregressGlobal','physioResidualVarRegressGlobal',...
    'designMatRespPulse','designMatResp','designMatPulse','rwdMeanRVT','rwdMeanRV',...
    'pulseBoldCorr');

plotColors = {[1 0 0], [0 0 1]};


i=0;

%% PULSE KERNEL
i=i+1; figure(i); clf
cols=2;
rows=ceil(length(roiNames)/cols);
for iRoi=1:length(roiNames)
    subplot(rows,cols,iRoi)
    for rwd=1:2
        hold on 
        plot(squeeze(mean(pulseKernel(:,iRoi,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
           for iSub=1:length(subFolders)
                pulseBoldCorr = corr(downsampledPulse{iSub,rwd},meanRoiTC{iSub,iRoi,rwd});
           end
    end
    title(['pulse kernel: ' roiNames{iRoi}]);
end

%%
i=i+1; figure(i); clf

for rwd=1:2
    hold on
    plot(squeeze(mean(rwdMeanPulseRate(:,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
end
title(['mean pulse']);


%%
i=i+1; figure(i); clf

for iRoi=1:length(roiNames)
    subplot(rows,cols,iRoi)
    for rwd=1:2
        for iSub=1:length(subFolders)
            meanRoiTrial(iSub,iRoi,rwd,:) = mean(reshape(meanRoiTC{iSub,iRoi,rwd},trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(meanRoiTrial(:,iRoi,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
    end
    title(['original trial: ' roiNames{iRoi}]);
end
%%
i=i+1; figure(i); clf

for iRoi=1:length(roiNames)
    subplot(rows,cols,iRoi)
    for rwd=1:2
        for iSub=1:length(subFolders)
            meanRoiTrialRegressPulse(iSub,iRoi,rwd,:) = mean(reshape(pulseResidualTC{iSub,iRoi,rwd},trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(meanRoiTrialRegressPulse(:,iRoi,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
    end
    title(['post pulse regression: ' roiNames{iRoi}]);
end

%%
i=i+1; figure(i); clf

for iRoi=1:length(roiNames)
    subplot(rows,cols,iRoi)
    for rwd=1:2
        for iSub=1:length(subFolders)
            pulseRegressedTC{iSub,rwd}(iRoi,:) = squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{iSub,rwd};
            pulseRegressedTrial(iSub,iRoi,rwd,:) = mean(reshape(pulseRegressedTC{iSub,rwd}(iRoi,:),trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(pulseRegressedTrial(:,iRoi,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
    end
    title(['regressed out: ' roiNames{iRoi}]);
end


%%
squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{iSub,rwd};