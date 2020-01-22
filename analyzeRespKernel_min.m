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
    'pulseBoldCorr','deconvLength',...
    'respKernel','respKernelRegress','respResidualTC','respResidualVar',...
    'respResidualTCregressGlobal','respResidualVarRegressGlobal',...
    'downsampledRV','downsampledRVT');

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
        plot(squeeze(mean(respKernel(:,iRoi,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
    end
    title(['resp kernel: ' roiNames{iRoi}]);
end

%%
i=i+1; figure(i); clf

for rwd=1:2
    hold on
    plot(squeeze(mean(rwdMeanRV(:,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
end
title(['mean RV']);


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
            meanRoiTrialRegressResp(iSub,iRoi,rwd,:) = mean(reshape(respResidualTC{iSub,iRoi,rwd},trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(meanRoiTrialRegressResp(:,iRoi,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
    end
    title(['post RV regression: ' roiNames{iRoi}]);
end

%%
i=i+1; figure(i); clf

for iRoi=1:length(roiNames)
    subplot(rows,cols,iRoi)
    for rwd=1:2
        for iSub=1:length(subFolders)
            respRegressedTC{iSub,rwd}(iRoi,:) = squeeze(respKernel(iSub,iRoi,rwd,:))'*designMatResp{iSub,rwd};
            respRegressedTrial(iSub,iRoi,rwd,:) = mean(reshape(respRegressedTC{iSub,rwd}(iRoi,:),trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(respRegressedTrial(:,iRoi,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
    end
    title(['regressed out: ' roiNames{iRoi}]);
end


%%
