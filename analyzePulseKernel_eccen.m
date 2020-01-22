close all
mrQuit
clear all
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';

load([dataFolder 'rwd_physioRegress_eccen.mat'], 'dataFolder', 'subFolders', 'numSubs','roiNames', ...
    'numRuns','numTRs','concatInfo',...
    'trialsPerRun','trialLength','junkedFrames','TR',...
    'ecgselect','ecgSampleRate','ecgTrial','ecgRunLength','ecgInterpMethod',...
    'numRuns','concatInfo','numTRs','globalMean',...
    'ecg','ecgPeaksDiff','ecgPeaksAmp','ecgRateTime','interpPeaksDiff',...
    'trialsPeakDiff','ecgPulseRate','interpPulseRate',...
    'trialsPulseRate','meanTrialPulseRateRun','meanTrialPeaksDiffRun',...
    'rwdMeanPeaksDiff','rwdMeanPulseRate','rwdPulseTC','downsampledPulse',...
    'rwdRvTC','rwdRvtTC','downsampledRV','downsampledRVT',...
    'designMatPulse',...
    'pulseKernel',...
    'pulseResidualTC',...
    'rv','rvt','respselect','resp','respPeaks','respPeaksDiff',...
    'respPeaksAmp','interpPeaksAmp','respRateTime',...
    'interpRespPeaksDiff','respTroughs','respTroughsAmp','interpTroughsAmp',...
    'physioKernel','physioResidualTC',...
    'designMatRespPulse','designMatResp','designMatPulse','rwdMeanRVT','rwdMeanRV',...
    'deconvLength',...
    'respKernel','respResidualTC',...
    'downsampledRV','downsampledRVT',...
    'eccen','ang','areas',...
    'numBinVoxels','binMeanTseries','subBinTrialResponse','binBorders','nbins');

plotColors = {[1 0 0], [0 0 1]};


i=0;

%% PULSE KERNEL
i=i+1; figure(i); clf
rows=1;
cols=nbins;%ceil(nbins/cols);
for ibin=1:nbins
    subplot(rows,cols,ibin)
    for rwd=1:2
%         for iSub=1:length(subFolders)
%             zscorePulseKernel(iSub,ibin,rwd,:) = zscore(squeeze(pulseKernel(iSub,ibin,rwd,:)));
%         end
        hold on 
%         plot(squeeze(mean(zscorePulseKernel(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
        plot(squeeze(mean(pulseKernel(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
        
%            for iSub=1:length(subFolders)
%                 pulseBoldCorr = corr(downsampledPulse{iSub,rwd},binMeanTseries{iSub,ibin,rwd});
%            end
    end
    
end
title(['pulse kernel']);
%%
i=i+1; figure(i); clf

for rwd=1:2
    hold on
    plot(squeeze(mean(rwdMeanPulseRate(:,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
end
title(['mean pulse']);


%%
i=i+1; figure(i); clf

for ibin=1:nbins
    subplot(rows,cols,ibin)
    for rwd=1:2
        for iSub=1:length(subFolders)
            meanRoiTrial(iSub,ibin,rwd,:) = mean(subBinTrialResponse{iSub,ibin,rwd},2);
        end
        hold on
        plot(squeeze(mean(meanRoiTrial(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
%         plot(,'color',plotColors{rwd},'linewidth', 3);
    end
    
end
title(['original trial']);
%%
i=i+1; figure(i); clf

for ibin=1:nbins
    subplot(rows,cols,ibin)
    for rwd=1:2
        for iSub=1:length(subFolders)
            meanBinTrialRegressPulse(iSub,ibin,rwd,:) = mean(reshape(pulseResidualTC{iSub,ibin,rwd},trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(meanBinTrialRegressPulse(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
    end
    
end
title(['post pulse regression']);
%%
i=i+1; figure(i); clf

for ibin=1:nbins
    subplot(rows,cols,ibin)
    for rwd=1:2
        for iSub=1:length(subFolders)
            pulseRegressedTC{iSub,rwd}(ibin,:) = squeeze(pulseKernel(iSub,ibin,rwd,:))'*designMatPulse{iSub,rwd};
            pulseRegressedTrial(iSub,ibin,rwd,:) = mean(reshape(pulseRegressedTC{iSub,rwd}(ibin,:),trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(pulseRegressedTrial(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
    end
    
end
title(['regressed out']);

%%
% squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{iSub,rwd};