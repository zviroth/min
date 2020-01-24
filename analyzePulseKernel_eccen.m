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
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;
linewidth=1;
markersize=10;
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
        plot(squeeze(mean(pulseKernel(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', linewidth);
        
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
    plot(squeeze(mean(rwdMeanPulseRate(:,rwd,:))),'color',plotColors{rwd},'linewidth', linewidth);
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
        plot(squeeze(mean(meanRoiTrial(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', linewidth);
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
        plot(squeeze(mean(pulseRegressedTrial(:,ibin,rwd,:))),'color',plotColors{rwd},'linewidth', linewidth);
    end
    
end
title(['regressed out']);

%%
i=i+1; figure(i); clf
% pulseKernel(iSub,ibin,rwd,timepoint)
pulseKernelStd = std(pulseKernel,0,4);
meanPulseKernelStd = squeeze(mean(pulseKernelStd));
stdPulseKernelStd = squeeze(std(pulseKernelStd));

for rwd=1:2
    dsErrorsurface(binCenters, squeeze(meanPulseKernelStd(:,rwd)), squeeze(stdPulseKernelStd(:,rwd))./sqrt(size(pulseKernelStd,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot(binCenters, squeeze(meanPulseKernelStd(:,rwd)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
end

    xlabel('eccentricity (deg)');
    ylabel('pulse kernel amplitude');
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
        'labelFontSize',7);
    axis square
    
%%
% squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{iSub,rwd};