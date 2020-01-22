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
goodSubs=2:length(subFolders);

i=0;

%% RESP KERNEL
i=i+1; figure(i); clf
cols=2;
rows=ceil(nbins/cols);
for ibin=1:nbins
    subplot(rows,cols,ibin)
    for rwd=1:2
        hold on 
        for iSub=1:length(subFolders)
            zscoreRespKernel(iSub,ibin,rwd,:) = zscore(squeeze(respKernel(iSub,ibin,rwd,:)));
        end
        plot(squeeze(mean(zscoreRespKernel(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
%            for iSub=1:length(subFolders)
%                 pulseBoldCorr = corr(downsampledPulse{iSub,rwd},binMeanTseries{iSub,ibin,rwd});
%            end
    end
    
end
title(['resp kernel']);
%%
i=i+1; figure(i); clf

for rwd=1:2
    hold on
    plot(squeeze(mean(rwdMeanRV(:,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
end
title(['mean RV']);


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
            meanBinTrialRegressResp(iSub,ibin,rwd,:) = mean(reshape(respResidualTC{iSub,ibin,rwd},trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(meanBinTrialRegressResp(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
    end
    
end
title(['post resp regression']);
%%
i=i+1; figure(i); clf

for ibin=1:nbins
    subplot(rows,cols,ibin)
    for rwd=1:2
        for iSub=1:length(subFolders)
            respRegressedTC{iSub,rwd}(ibin,:) = squeeze(respKernel(iSub,ibin,rwd,:))'*designMatResp{iSub,rwd};
            respRegressedTrial(iSub,ibin,rwd,:) = mean(reshape(respRegressedTC{iSub,rwd}(ibin,:),trialLength,[]),2);
        end
        hold on
        plot(squeeze(mean(respRegressedTrial(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', 3);
    end
    
end
title(['regressed out']);

%%
% squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{iSub,rwd};