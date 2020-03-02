close all
mrQuit
clear all
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
deconvLength=12;
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
    'respResidualTCregressGlobal','respResidualVarRegressGlobal');

plotColors = {[1 0 0], [0 0 1]};
plotStyles = {'-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.'};
linewidth=2;
i=0;
nkernels = size(physioKernel,4)/deconvLength;%RV, RVT, pulse
kernelStr = {'RV','RVT','pulse'};
%%
i=i+1; figure(i); clf
cols=1;
rows=3;
subplot(rows,cols,1) %RV
for rwd=1:2
%     plot(squeeze(rwdMeanRV(:,rwd,:))','color',plotColors{rwd},'linewidth', linewidth-1);
    hold on
    plot(squeeze(mean(rwdMeanRV(:,rwd,:))),'color',plotColors{rwd},'linewidth', linewidth);
end
title(['mean RV']);
subplot(rows,cols,2) %RVT
for rwd=1:2
%     plot(squeeze(rwdMeanRVT(:,rwd,:))','color',plotColors{rwd},'linewidth', linewidth-1);
    hold on
    plot(squeeze(mean(rwdMeanRVT(:,rwd,:))),'color',plotColors{rwd},'linewidth', linewidth);
end
title(['mean RVT']);
subplot(rows,cols,3) %pulse
for rwd=1:2
    hold on
    plot(squeeze(mean(rwdMeanPulseRate(:,rwd,:))),'color',plotColors{rwd},'linewidth', linewidth);
end
title(['mean pulse']);

%% PHYSIO KERNEL
i=i+1; figure(i); clf
cols=length(roiNames);
rows=nkernels;

for ikernel=1:nkernels
    for iRoi=1:length(roiNames)
        subplot(rows,cols,iRoi + cols*(ikernel-1))
        for rwd=1:2
            hold on
            plot(squeeze(mean(physioKernel(:,iRoi,rwd,(ikernel-1)*deconvLength+1:ikernel*deconvLength))),'color',plotColors{rwd},'linewidth', linewidth);
%             for iSub=1:length(subFolders)
%                 pulseBoldCorr = corr(downsampledPulse{iSub,rwd},meanRoiTC{iSub,iRoi,rwd});
%             end
        end
    end
    title([kernelStr{ikernel} ' kernel: ' roiNames{iRoi}]);
end

%%
i=i+1; figure(i); clf
rows=2;
cols = ceil(length(subFolders)/rows);
iRoi=2;
ikernel=2;
for iSub=1:length(subFolders)
    subplot(rows,cols,iSub)
   for rwd=1:2
       plot(squeeze(physioKernel(iSub,iRoi,rwd,(ikernel-1)*deconvLength+1:ikernel*deconvLength)),'color',plotColors{rwd},'linewidth', linewidth);
       hold on;
   end
end


%%
i=i+1; figure(i); clf
rows=1;
for iRoi=1:length(roiNames)
    subplot(rows,cols,iRoi)
    for rwd=1:2
        for iSub=1:length(subFolders)
            meanRoiTrial(iSub,iRoi,rwd,:) = mean(reshape(meanRoiTC{iSub,iRoi,rwd},trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(meanRoiTrial(:,iRoi,rwd,:))),'color',plotColors{rwd},'linewidth', linewidth);
    end
    title(['original trial: ' roiNames{iRoi}]);
end
%%
i=i+1; figure(i); clf
rows=1;

for iRoi=1:length(roiNames)
    %     subplot(rows,cols,iRoi)
    subplot(rows,cols,iRoi)
    for rwd=1:2
        for iSub=1:length(subFolders)
            meanRoiTrialRegressPhysio(iSub,iRoi,rwd,:) = mean(reshape(physioResidualTC{iSub,iRoi,rwd},trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(meanRoiTrialRegressPhysio(:,iRoi,rwd,:))),'color',plotColors{rwd},'linewidth', linewidth);
    end
    
    title(['post regression: ' roiNames{iRoi}]);
end

%%
i=i+1; figure(i); clf
rows=nkernels+1;

for iRoi=1:length(roiNames)
    for rwd=1:2
        for ikernel=1:nkernels
            subplot(rows,cols,iRoi + (ikernel-1)*cols)
            for iSub=1:length(subFolders)
                %                 designMatRespPulse
                physioRegressedTC{iSub,rwd}(iRoi,ikernel,:) = squeeze(physioKernel(iSub,iRoi,rwd,(ikernel-1)*deconvLength+1:ikernel*deconvLength))'*designMatRespPulse{iSub,rwd}((ikernel-1)*deconvLength+1:ikernel*deconvLength,:);
                %                 pulseRegressedTC{iSub,rwd}(iRoi,:) = squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{iSub,rwd};
                %                 pulseRegressedTrial(iSub,iRoi,rwd,:) = mean(reshape(pulseRegressedTC{iSub,rwd}(iRoi,:),trialLength,[]),2);
                physioRegressedTrial(iSub,iRoi,rwd,ikernel,:) = mean(reshape(physioRegressedTC{iSub,rwd}(iRoi,ikernel,:),trialLength,[]),2);
            end
            hold on
            plot(squeeze(mean(physioRegressedTrial(:,iRoi,rwd,ikernel,:))),'color',plotColors{rwd},'linewidth', linewidth);
        end
        for iSub=1:length(subFolders)
            physioRegressedTC{iSub,rwd}(iRoi,nkernels+1,:) = squeeze(physioKernel(iSub,iRoi,rwd,:))'*designMatRespPulse{iSub,rwd};
            physioRegressedTrial(iSub,iRoi,rwd,nkernels+1,:) = mean(reshape(physioRegressedTC{iSub,rwd}(iRoi,nkernels+1,:),trialLength,[]),2);
        end
        subplot(rows,cols,iRoi + (nkernels)*cols); hold on
        plot(squeeze(mean(physioRegressedTrial(:,iRoi,rwd,nkernels+1,:))),'color',plotColors{rwd},'linewidth', linewidth+1);
    end
    
    title(['regressed out: ' roiNames{iRoi}]);
end


%%
squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{iSub,rwd};