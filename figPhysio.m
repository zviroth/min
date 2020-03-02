close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 0;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';

onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
elseif onlyCorrect==0
    onlyCorrectString = '_validresponse';
end
zScoreString = '';
fmriUnits = '% change image intensity';
if toZscore
    zScoreString = '_zscored';
    fmriUnits = 'std image intensity';
end

globalMeanString = '';
if regressGlobalMean
    globalMeanString = '_globalRegressed';
end
ConcatProjStr = '';
if ConcatProj
    ConcatProjStr = 'ConcatProj';
end
load([dataFolder 'rwdTC_physio' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'globalMean','regressBetasGlobal','runRwd',...
    'subRoiRuns','runMeanFFT',...
    'allVoxTrialResponse','allVoxTaskPhase','allVoxTaskAmp','allVoxTaskCo',...
    'voxTrials','voxGoodTrials','meanVoxTrial',...
    'maxRT',...
    'ecgselect','ecgSampleRate','ecgTrial','ecgRunLength','ecgInterpMethod',...
    'ecg','ecgPulseRate','interpPulseRate',...
    'respselect','resp',...
    'rwdPulseTC','rwdRvTC',...
    'designMatPulse','designMatRespPulse','designMatResp','deconvLength',...
    'allGoodTrials','allGoodTRs');

plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
lineLength = 0.2;
lineWidth = 2;
markerSize = 10;
fontsize=9;
dsSurfaceContrast = 1;
dsSurfaceAlpha = 0.15;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
% ROIs = [2];
eccMin = 0.2;
eccMax = 70;
nbins = 12;
TR=1.5;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end

goodSubs = [1:length(subFolders)]; 

%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;

%% mean + STD of thresholds
subThresh = mean(subMeanThresh,2);%mean over rwd, which is already averaged over staircases and runs
mean(subThresh);
std(subThresh);


%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd

for iSub = 1:length(goodSubs)%1:length(goodSubs)%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(roiNames)
        for rwd=1:2
%             reshapedTrials = reshape(subTrialResponse{goodSubs(iSub),iRoi,rwd},10,[]);



        trRV{iSub,rwd} = reshape(rwdRvTC{goodSubs(iSub),rwd}, ecgTrial,[]);%this is only good trials!
        subMeanRV(iSub,rwd,:) = nanmean(trRV{iSub,rwd},2);
        subRespStd(iSub,rwd) = std(subMeanRV(iSub,rwd,:));%std amplitude of mean
        subRespVar(iSub,rwd,:) = nanstd(trRV{iSub,rwd},0,2);%timepoint variability
        respStd = nanstd(trRV{iSub,rwd});%std amp per trial
        subRespStdVar(iSub,rwd) = std(respStd);
        f = fft(trRV{iSub,rwd});
        subRespPhVar(iSub,rwd) = nanstd(angle(f(2,:)));
        subRespAmpVar(iSub,rwd) = nanstd(abs(f(2,:)));
        
        
        trPulse{iSub,rwd} = reshape(rwdPulseTC{goodSubs(iSub),rwd}, ecgTrial,[]);%this is only good trials!
        subMeanPulse(iSub,rwd,:) = nanmean(trPulse{iSub,rwd},2);
        subPulseStd(iSub,rwd) = std(subMeanPulse(iSub,rwd,:));%std amplitude of mean
        subPulseVar(iSub,rwd,:) = nanstd(trPulse{iSub,rwd},0,2);%timepoint variability
        pulseStd = nanstd(trPulse{iSub,rwd});%std amp per trial
        subPulseStdVar(iSub,rwd) = std(pulseStd);
        f = fft(trPulse{iSub,rwd});
        subPulsePhVar(iSub,rwd) = nanstd(angle(f(2,:)));
        subPulseAmpVar(iSub,rwd) = nanstd(abs(f(2,:)));
        
            %REGRESS OUT PHYSIO
            %roiMeanTseries{iSub,iRoi,rwd} includes all trials
           designMatPulse{goodSubs(iSub),rwd}= designMatPulse{goodSubs(iSub),rwd}*ecgSampleRate;%change to beats per sec
           pulseKernel(iSub,iRoi,rwd,:) = designMatPulse{goodSubs(iSub),rwd}'\subTrialResponse{goodSubs(iSub),iRoi,rwd}(:);%this is only good trials!
           pulseResidualTC{iSub,iRoi,rwd} = subTrialResponse{goodSubs(iSub),iRoi,rwd}(:)' - squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{goodSubs(iSub),rwd};

           designMatResp{goodSubs(iSub),rwd} = designMatResp{goodSubs(iSub),rwd}*ecgSampleRate;%change to beats per sec
           respKernel(iSub,iRoi,rwd,:) = designMatResp{goodSubs(iSub),rwd}'\subTrialResponse{goodSubs(iSub),iRoi,rwd}(:);
           respResidualTC{iSub,iRoi,rwd} = subTrialResponse{goodSubs(iSub),iRoi,rwd}(:)' - squeeze(respKernel(iSub,iRoi,rwd,:))'*designMatResp{goodSubs(iSub),rwd};

           designMatRespPulse{goodSubs(iSub),rwd} = designMatRespPulse{goodSubs(iSub),rwd}*ecgSampleRate;%change to beats per sec
           physioKernel(iSub,iRoi,rwd,:) = designMatRespPulse{goodSubs(iSub),rwd}'\subTrialResponse{goodSubs(iSub),iRoi,rwd}(:);
           physioResidualTC{iSub,iRoi,rwd} = subTrialResponse{goodSubs(iSub),iRoi,rwd}(:)' - squeeze(physioKernel(iSub,iRoi,rwd,:))'*designMatRespPulse{goodSubs(iSub),rwd};
           subTrialResponse{goodSubs(iSub),iRoi,rwd} = reshape(physioResidualTC{iSub,iRoi,rwd}(:), trialLength, length(physioResidualTC{iSub,iRoi,rwd}(:))/trialLength);


        end
     end
end

% regress physio per bin
iRoi=2;
for iSub = 1:length(goodSubs)%1:length(goodSubs)%length(subdirs)
    for rwd=1:2
        for ibin=1:nbins
            binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
            numBinVoxels(iSub,ibin) = sum(binVoxels);
            %                 binVoxels = binVoxels & areas{iSub,iRoi}==1;%ONLY V1
            binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
            binMeanTseries{iSub,ibin,rwd} = nanmean(roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,allGoodTRs{iSub,rwd}));%mean timecourse across voxels
            subBinTrialResponse{iSub,ibin,rwd} = reshape(binMeanTseries{iSub,ibin,rwd}, trialLength, length(binMeanTseries{goodSubs(iSub),ibin,rwd})/trialLength);
            
            
            binPulseKernel(iSub,ibin,rwd,:) = designMatPulse{goodSubs(iSub),rwd}'\binMeanTseries{iSub,ibin,rwd}';
            binPulseResidualTC{iSub,ibin,rwd} = binMeanTseries{iSub,ibin,rwd} - squeeze(binPulseKernel(iSub,ibin,rwd,:))'*designMatPulse{goodSubs(iSub),rwd};
            
            binPhysioKernel(iSub,ibin,rwd,:) = designMatRespPulse{goodSubs(iSub),rwd}'\binMeanTseries{iSub,ibin,rwd}';
            binPhysioResidualTC{iSub,ibin,rwd} = binMeanTseries{iSub,ibin,rwd} - squeeze(binPhysioKernel(iSub,ibin,rwd,:))'*designMatRespPulse{goodSubs(iSub),rwd};
            
            binRespKernel(iSub,ibin,rwd,:) = designMatResp{goodSubs(iSub),rwd}'\binMeanTseries{iSub,ibin,rwd}';
            binRespResidualTC{iSub,ibin,rwd} = binMeanTseries{iSub,ibin,rwd} - squeeze(binRespKernel(iSub,ibin,rwd,:))'*designMatResp{goodSubs(iSub),rwd};
        end 
    end
end


for iSub = 1:length(goodSubs)%length(subdirs)
    for rwd=1:2
        temp = trialCorrectness{goodSubs(iSub),rwd}(:,2:end-1);%exclude the first trial and the last trial
        trialCorrectVec{iSub,rwd} = temp(:); %include trials with no response
        subPropCorrect(iSub,rwd) = sum(trialCorrectVec{iSub,rwd})/length(trialCorrectVec{iSub,rwd});
        
        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!
        temp = trialResponse{goodSubs(iSub),rwd}(:,2:end-1);
%         temp = trialResponse{goodSubs(iSub),rwd}(:);
        trialResponseVec = temp(:);
%         temp = trialRT{goodSubs(iSub),rwd}(:);
        temp = trialRT{goodSubs(iSub),rwd}(:,2:end-1);
        allTrialRTvec = temp(:);
        
        if onlyCorrect==0
            trialRTvec{iSub,rwd} = temp(trialResponseVec>0 & allTrialRTvec>0 & allTrialRTvec<maxRT);%only trials with a response, we're assuming onlyCorrect==0
        elseif onlyCorrect==1 %only correct
            trialRTvec{iSub,rwd} = temp(trialCorrectVec{iSub,rwd}==1 & allTrialRTvec>0 & allTrialRTvec<maxRT);
        elseif onlyCorrect==2 %only incorrect
            trialRTvec{iSub,rwd} = temp(trialCorrectVec{iSub,rwd}==0 & trialResponseVec>0 & allTrialRTvec>0 & allTrialRTvec<maxRT);
        else
            trialRTvec{iSub,rwd} = temp(:);
        end

    end

    %pulse & respiration deconvolution design matrices
    temp = [designMatPulse{goodSubs(iSub),1} designMatPulse{goodSubs(iSub),2}];
    designMatPulseTrials = reshape(temp,deconvLength,trialLength,[]);%10 deconvolution points, 10 trial timepoints, X trials
    temp = [designMatResp{goodSubs(iSub),1} designMatResp{goodSubs(iSub),2}];
    designMatRespTrials = reshape(temp,deconvLength,trialLength,[]);%10 deconvolution points, 10 trial timepoints, X trials
    temp = [designMatRespPulse{goodSubs(iSub),1} designMatRespPulse{goodSubs(iSub),2}];
    designMatPhysioTrials = reshape(temp,deconvLength,trialLength,[]);%10 deconvolution points, 10 trial timepoints, X trials

end

realSubPulseKernelStd = std(pulseKernel,0,4);
realSubRespKernelStd = std(respKernel,0,4);


%%
load(fullfile(dataFolder,'subColors.mat'),'subColor', 'scatterCmap', 'colorSubFolders', 'colorGoodSubs');
linewidth = 1;
markersize=10;
subjects = size(subPropCorrect,1);
rewards = size(subPropCorrect,2);
[rwdNum, subNum] = meshgrid(1:rewards, 1:subjects);
%     scatterCmap = scatterCmap(end:-1:1,:);%invert color map
colormap(scatterCmap);


rows=2;
cols = 8;
subplots = {1:3, 4:6, 7 , 9, cols+1:cols+3, cols+4:cols+6, cols+7 , cols+9};
% subplots = {1 3 5};
lineLength = 0.2;
lineWidth = 2;
markerSize = 10;

%% RESPIRATION KERNEL
i=0;

i=i+1; figure(i); clf

% pulseKernel(iSub,ibin,rwd,timepoint)
respKernelStd = std(binRespKernel,0,4);
meanRespKernelStd = squeeze(mean(respKernelStd));
stdRespKernelStd = squeeze(std(respKernelStd));
respKernelStdDiff = squeeze(respKernelStd(:,:,1) - respKernelStd(:,:,2));
meanRespKernelStdDiff = squeeze(mean(respKernelStdDiff));
stdRespKernelStdDiff = squeeze(std(respKernelStdDiff));

subplot(rows,cols,subplots{1})
iRoi=2;
meanRespKernel = squeeze(mean(respKernel(:,iRoi,:,:)));
stdRespKernel = squeeze(std(respKernel(:,iRoi,:,:)));
for rwd=1:2
    dsErrorsurface(TR*(1:trialLength), squeeze(meanRespKernel(rwd,:)), squeeze(stdRespKernel(rwd,:))./sqrt(size(respKernel,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot(TR*(1:trialLength), squeeze(meanRespKernel(rwd,:)),'-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
end
% ylabel('RV kernel');
t = {'BOLD signal'; ['(' fmriUnits ')'] };
ylabel(t);

subplot(rows,cols,subplots{2})
for rwd=1:2
    dsErrorsurface(binCenters, squeeze(meanRespKernelStd(:,rwd)), squeeze(stdRespKernelStd(:,rwd))./sqrt(size(respKernelStd,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot(binCenters, squeeze(meanRespKernelStd(:,rwd)),'.','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
end
% ylabel('respiration kernel amplitude');
t = {'respiration kernel amplitude'; ['(' fmriUnits ')'] };
ylabel(t);


% subplot(rows,cols,subplots{2})
%     dsErrorsurface(binCenters, meanRespKernelStdDiff, stdRespKernelStdDiff./sqrt(size(respKernelStd,1)), [0 0 0],dsSurfaceAlpha);
%     hold on
%     plot(binCenters, meanRespKernelStdDiff,'.-','Color', [0 0 0],'linewidth',linewidth,'markersize',markersize);
% ylabel('\Delta respiration kernel amplitude');

for isubplot=1:2
    subplot(rows,cols,subplots{isubplot})
    xlabel('eccentricity (deg)');
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -10/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
    'labelFontSize',fontsize);
    axis square
end


% subject bar plot - respiration kernel amplitude
iRoi=2;
subplot(rows,cols,subplots{3})
% Resp kernel std
for iSub=1:length(subFolders)
    %find the correct color
    for j=1:length(colorGoodSubs)
        if strcmp(subFolders{iSub},colorSubFolders{colorGoodSubs(j)})
           newSubColor(iSub,:) = subColor(j,:);
        end
    end
%     subColor(iSub,:) = scatterCmap(1+floor((smallSubDiff(iSub) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
    plot(1:rewards,squeeze(realSubRespKernelStd(iSub,iRoi,:)),'Color',newSubColor(iSub,:),'linewidth',linewidth);
    hold on
end
temp = realSubRespKernelStd(:,iRoi,:);
scatter(rwdNum(:),temp(:),markerSize, [newSubColor; newSubColor] ,'filled');
% MEAN
for irwd=1:rewards
    scatter(irwd,squeeze(mean(realSubRespKernelStd(:,iRoi,irwd))),markerSize*4,[0 0 0],'filled');
    plot(1:rewards,squeeze(mean(realSubRespKernelStd(:,iRoi,:))),'Color',[0 0 0],'linewidth', 2*linewidth);
end


%
for isubplot=3
    subplot(rows,cols,subplots{isubplot})
    switch isubplot
        case 3
            ylabel('respiration kernel amplitude');

    end 
    xlabel('reward');
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -38/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
        'xTick',[1 2], 'xTickLabel', {'high','low'},...
        'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',fontsize,...
        'yAxisMajorTickLen',-4/32);
end
set(gcf,'position',[10 10 24 	14]);
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig_respKernel_' ConcatProjStr '.pdf']);


%% RESPIRATION VOLUME

for iSub=1:length(subFolders)
    for rwd=1:2
       temp = reshape(rwdRvTC{iSub,rwd},ecgTrial,[]);
       subRespTrial(iSub,rwd,:) = nanmean(temp,2);
    end
end
subRespTrial = subRespTrial*ecgSampleRate;%changing to beats-per-second
meanRespTrial = squeeze(nanmean(subRespTrial));
stdRespTrial = squeeze(nanstd(subRespTrial));
subRespDiff = squeeze(subRespTrial(:,1,:) - subRespTrial(:,2,:));
meanRespDiff = squeeze(nanmean(subRespDiff));
stdRespDiff = squeeze(nanstd(stdRespTrial));


i=i+1; figure(i); clf
%mean pulse kernel
subplot(rows,cols,subplots{1})
for rwd=1:2
    dsErrorsurface((0:ecgTrial-1)/ecgSampleRate, meanRespTrial(rwd,:),stdRespTrial(rwd,:)./sqrt(size(rwdPulseTC,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot((0:ecgTrial-1)/ecgSampleRate,meanRespTrial(rwd,:), 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',20);
end
ylabel('RV');

%difference in pulse rate between runs
subplot(rows,cols,subplots{2})

dsErrorsurface((0:ecgTrial-1)/ecgSampleRate, meanRespDiff,stdRespDiff./sqrt(size(rwdPulseTC,1)), [0  0 0],dsSurfaceAlpha);
hold on
plot((0:ecgTrial-1)/ecgSampleRate,meanRespDiff, 'Color', [0 0 0], 'linewidth', linewidth,'markersize',20);
hline(0)
ylabel('\Delta RV');

for isubplot=1:2
    subplot(rows,cols,subplots{isubplot})
    xlabel('time (s)');
        drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -10/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
    'labelFontSize',fontsize);
    axis square
end


% scatter plot of subjects respiration
subplot(rows,cols,subplots{3})

for iSub=1:length(subFolders)
%     subColor(iSub,:) = scatterCmap(1+floor((smallSubDiff(iSub) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
    plot(1:rewards,subRespStd(iSub,:),'Color',newSubColor(iSub,:),'linewidth',linewidth);
    hold on
end
scatter(rwdNum(:),subRespStd(:),markerSize, [newSubColor; newSubColor] ,'filled');
% MEAN
for irwd=1:rewards
    scatter(irwd,mean(subRespStd(:,irwd)),markerSize*4,[0 0 0],'filled');
    plot(1:rewards,mean(subRespStd),'Color',[0 0 0],'linewidth', 2*linewidth);
end


xlabel('reward');
ylabel('RV amplitude ');
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -50/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
        'xTick',[1 2], 'xTickLabel', {'high','low'},...
        'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',fontsize,...
        'yAxisMajorTickLen',-4/32);
set(gcf,'position',[10 10 24 	14]);

print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig_respRate.pdf']);



%% PULSE KERNEL




% pulseKernel(iSub,ibin,rwd,timepoint)
pulseKernelStd = std(binPulseKernel,0,4);
meanPulseKernelStd = squeeze(mean(pulseKernelStd));
stdPulseKernelStd = squeeze(std(pulseKernelStd));
pulseKernelStdDiff = squeeze(pulseKernelStd(:,:,1) - pulseKernelStd(:,:,2));
meanPulseKernelStdDiff = squeeze(mean(pulseKernelStdDiff));
stdPulseKernelStdDiff = squeeze(std(pulseKernelStdDiff));



i=i+1; figure(i); clf
%mean pulse kernel
subplot(rows,cols,subplots{1})
iRoi=2;
meanPulseKernel = squeeze(mean(pulseKernel(:,iRoi,:,:)));
stdPulseKernel = squeeze(std(pulseKernel(:,iRoi,:,:)));
for rwd=1:2
    dsErrorsurface(TR*(1:trialLength), squeeze(meanPulseKernel(rwd,:)), squeeze(stdPulseKernel(rwd,:))./sqrt(size(pulseKernel,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot(TR*(1:trialLength), squeeze(meanPulseKernel(rwd,:)),'-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
end
t = {'BOLD signal'; ['(' fmriUnits ')'] };
ylabel(t);
ylim([-.02 0.03]);
% ylabel('pulse kernel');
xlabel('time (sec)');
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -14/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
    'labelFontSize',fontsize);

%pulse kernel amplitude vs eccentricity
subplot(rows,cols,subplots{2})
for rwd=1:2
    dsErrorsurface(binCenters, squeeze(meanPulseKernelStd(:,rwd)), squeeze(stdPulseKernelStd(:,rwd))./sqrt(size(pulseKernelStd,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot(binCenters, squeeze(meanPulseKernelStd(:,rwd)),'.','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
end
t = {'pulse kernel amplitude'; ['(' fmriUnits ')'] };
ylabel(t);
% ylabel('pulse kernel amplitude');
xlabel('eccentricity (deg)');
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -8/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
    'labelFontSize',fontsize);

% subplot(rows,cols,subplots{2})
%     dsErrorsurface(binCenters, meanPulseKernelStdDiff, stdPulseKernelStdDiff./sqrt(size(pulseKernelStd,1)), [0 0 0],dsSurfaceAlpha);
%     hold on
%     plot(binCenters, meanPulseKernelStdDiff,'.-','Color', [0 0 0],'linewidth',linewidth,'markersize',markersize);
% ylabel('\Delta pulse kernel amplitude');

for isubplot=1:2
    subplot(rows,cols,subplots{isubplot})
%     drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
%     'labelFontSize',fontsize);
    axis square
end

% subject bar plot - pulse kernel amplitude
iRoi=2;
subplot(rows,cols,subplots{3})
% Pulse kernel std
l = size(scatterCmap,1);
for iSub=1:length(subFolders)
    %find the correct color
    for j=1:length(colorGoodSubs)
        if strcmp(subFolders{iSub},colorSubFolders{colorGoodSubs(j)})
           newSubColor(iSub,:) = subColor(j,:);
        end
    end
%     subColor(iSub,:) = scatterCmap(1+floor((smallSubDiff(iSub) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
    plot(1:rewards,squeeze(realSubPulseKernelStd(iSub,iRoi,:)),'Color',newSubColor(iSub,:),'linewidth',linewidth);
    hold on
end
temp = realSubPulseKernelStd(:,iRoi,:);
scatter(rwdNum(:),temp(:),markerSize, [newSubColor; newSubColor] ,'filled');
% MEAN
for irwd=1:rewards
    scatter(irwd,squeeze(mean(realSubPulseKernelStd(:,iRoi,irwd))),markerSize*4,[0 0 0],'filled');
    plot(1:rewards,squeeze(mean(realSubPulseKernelStd(:,iRoi,:))),'Color',[0 0 0],'linewidth', 2*linewidth);
end


%
for isubplot=3
    subplot(rows,cols,subplots{isubplot})
    switch isubplot
        case 3
            ylabel('pulse kernel amplitude');

    end 
    xlabel('reward');
    if ConcatProj
        drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -26/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
        'xTick',[1 2], 'xTickLabel', {'high','low'},...
        'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',fontsize,...
        'yAxisMajorTickLen',-4/32);
    else
        drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -62/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
        'xTick',[1 2], 'xTickLabel', {'high','low'},...
        'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',fontsize,...
        'yAxisMajorTickLen',-4/32);
    end
end
set(gcf,'position',[10 10 24 	14]);
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig_pulseKernel_' ConcatProjStr '.pdf']);

%%
%% PULSE RATE

for iSub=1:length(subFolders)
    for rwd=1:2
       temp = reshape(rwdPulseTC{iSub,rwd},ecgTrial,[]);
       subPulseTrial(iSub,rwd,:) = nanmean(temp,2);
    end
end
subPulseTrial = subPulseTrial*ecgSampleRate;%changing to beats-per-second
meanPulseTrial = squeeze(nanmean(subPulseTrial));
stdPulseTrial = squeeze(nanstd(subPulseTrial));
subPulseDiff = squeeze(subPulseTrial(:,1,:) - subPulseTrial(:,2,:));
meanPulseDiff = squeeze(nanmean(subPulseDiff));
stdPulseDiff = squeeze(nanstd(subPulseDiff));

i=i+1; figure(i); clf

%mean pulse rate 
subplot(rows,cols,subplots{1})
for rwd=1:2
    dsErrorsurface((0:ecgTrial-1)/ecgSampleRate, meanPulseTrial(rwd,:),stdPulseTrial(rwd,:)./sqrt(size(rwdPulseTC,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot((0:ecgTrial-1)/ecgSampleRate,meanPulseTrial(rwd,:), 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',20);
end
ylabel('pulse (beats/sec)');

%difference in pulse rate between runs
subplot(rows,cols,subplots{2})

dsErrorsurface((0:ecgTrial-1)/ecgSampleRate, meanPulseDiff,stdPulseDiff./sqrt(size(rwdPulseTC,1)), [0  0 0],dsSurfaceAlpha);
hold on
plot((0:ecgTrial-1)/ecgSampleRate,meanPulseDiff, 'Color', [0 0 0], 'linewidth', linewidth,'markersize',20);
hline(0);
ylabel('\Delta pulse (beats/sec)');

for isubplot=1:2
    subplot(rows,cols,subplots{isubplot})
    xlabel('time (s)');
        drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
    'labelFontSize',fontsize);
    axis square
end


% scatter plot of subjects pulse
subplot(rows,cols,subplots{3})

for iSub=1:length(subFolders)
%     subColor(iSub,:) = scatterCmap(1+floor((smallSubDiff(iSub) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
    plot(1:rewards,subPulseStd(iSub,:),'Color',newSubColor(iSub,:),'linewidth',linewidth);
    hold on
end
scatter(rwdNum(:),subPulseStd(:),markerSize, [newSubColor; newSubColor] ,'filled');
% MEAN
for irwd=1:rewards
    scatter(irwd,mean(subPulseStd(:,irwd)),markerSize*4,[0 0 0],'filled');
    plot(1:rewards,mean(subPulseStd),'Color',[0 0 0],'linewidth', 2*linewidth);
end


xlabel('reward');
ylabel('pulse amplitude ');
drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -54/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
    'xTick',[1 2], 'xTickLabel', {'high','low'},...
    'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',fontsize,...
    'yAxisMajorTickLen',-4/32);
set(gcf,'position',[10 10 24 	14]);

print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig_pulseRate.pdf']);

