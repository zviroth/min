close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
regressGlobalMean = 0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
TR=1.5;

onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
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
load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'globalMean','regressBetasGlobal','runRwd');

% trialLength=10;
clear subMeanResponse trialStd meanTrialStd runStd meanRunStd subRwdStd trialFFTamp trialFFTphase meanTrialFFTamp meanTrialFFTphase 
clear stdTrialFFTphase runFFTamp runFFTphase meanRunFFTamp meanRunFFTphase stdRunFFTphase subFFTamp subFFTphase
clear trialMinMax meanTrialMinMax runMinMax meanRunMinMax subRwdMinMax legendCellArray
clear stdRunMax stdTrialMax
clear groupLabels pVal
clear yError yMin yMax permDiff meanPermDiff
clear permSubDiff realSubDiff permSubResp realSubResp meanPermTC permSubMeanTC realSubMeanTC meanRealTC
plotColors = { [1 0 0],[0 0 1],[0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
ROIs = [1 2 7 8];
ROIs = [2];

dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;



%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;



%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd binCenters


goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% subjects. 0.0457  0.0140,     0.0553    0.0214

for iSub = goodSubs%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(roiNames)
        for rwd=1:2
            reshapedTrials = reshape(subTrialResponse{iSub,iRoi,rwd},10,[]);
            if iSub==goodSubs(1)
                allTrials{iRoi,rwd} = reshapedTrials;
            else
                allTrials{iRoi,rwd} = [allTrials{iRoi,rwd} reshapedTrials];
            end
        end
     end
end





%% bin voxels within each ROI according to eccentricity, then average within bins
%bins by log eccentricity
eccMin = 0.2;
eccMax = 70;
nbins = 12;
% binBorders = (linspace(eccMin^3,eccMax^3,nbins+1))^(1/3)
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);

nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end

clear subBinTrialResponse subBinResponse binVoxels numBinVoxels meanSubVoxCorr meanVoxCorr
for iSub = 1:length(goodSubs)%length(subdirs)
    for iRoi= ROIs
        for rwd=1:2
            temp = trialCorrectness{goodSubs(iSub),rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);
            temp = trialResponse{goodSubs(iSub),rwd}(:,2:end-1);
            trialResponseVec = temp(:);
            if onlyCorrect ==1 %ONLY CORRECT
                    goodTrials = trialCorrectnessVec==1;
                    %                 numTrials = sum(trialCorrectnessVec);
                elseif onlyCorrect ==2 % ONLY INCORRECT!!!
                    goodTrials = trialCorrectnessVec==0 & trialResponseVec>0;
                else % including all trials with a response
                    goodTrials = trialResponseVec>0;
            end
                
            for ibin=1:nbins
                binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
%                 binVoxels = binVoxels & areas{iSub,iRoi}==1;%ONLY V1
                binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries = nanmean(roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:));%mean timecourse across voxels
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = reshape(binMeanTseries, trialLength, length(binMeanTseries)/trialLength);
                
                
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = subBinTrialResponse{iSub,iRoi,ibin,rwd}(:,goodTrials);
                
                
                reshapedTrials = reshape(subBinTrialResponse{iSub,iRoi,ibin,rwd},trialLength,[]);
                if iSub==1
                    allBinTrials{iRoi,ibin,rwd} = reshapedTrials;
                else
                    allBinTrials{iRoi,ibin,rwd} = [allBinTrials{iRoi,ibin,rwd} reshapedTrials];
                end
                %average per subject
                subBinResponse(iSub,iRoi,ibin,rwd,:) = mean(subBinTrialResponse{iSub,iRoi,ibin,rwd},2);
                %trial-by-trial variability per subject
                subBinVar(iSub,iRoi,ibin,rwd) = mean(std(subBinTrialResponse{iSub,iRoi,ibin,rwd},0,2));
                
                temp = fft(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                roiBinFftAmp = abs(temp(2,:));
                roiBinFftPh = angle(temp(2,:));
                subBinFftAmpVar(iSub,iRoi,ibin,rwd) = std(roiBinFftAmp);
                subBinFftPhVar(iSub,iRoi,ibin,rwd) = circ_std(roiBinFftPh');
                
                singleTrialStd = std(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                subBinStdVar(iSub,iRoi,ibin,rwd) = std(singleTrialStd);
                subBinStdMean(iSub,iRoi,ibin,rwd) = mean(singleTrialStd);

            end
                %covariation between voxels
                temp = repmat(goodTrials', trialLength,1);
                goodTimepoints = temp(:);
%                 goodTC = roiTC{iSub,iRoi,rwd}.tSeries(binVoxels,goodTimepoints);
                goodTC = roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(:,goodTimepoints);
                %   maybe remove mean trial per voxel?
                temp=triu(corr(goodTC'));
                meanSubVoxCorr(iSub,iRoi,rwd) = nanmean(temp(:));
        end
    end
end
meanVoxCorr = squeeze(nanmean(meanSubVoxCorr));
fisherZ = atan(meanSubVoxCorr);
for iRoi= ROIs
    [h,  roiCorrPval(iRoi)] = ttest(fisherZ(:,iRoi,1)-fisherZ(:,iRoi,2));

end

subBinAmp = std(subBinResponse,0,5);%get mean timecourse amplitude per subject.(iRoi,ibin,rwd,1)
binMeanAmp = squeeze(mean(subBinAmp,1)); %mean over subjects. (iRoi,ibin,rwd)
binStdAmp = squeeze(std(subBinAmp,0,1)); %(iRoi,ibin,rwd)
binDiffAmp = squeeze(binMeanAmp(:,:,1) - binMeanAmp(:,:,2));%(iRoi,ibin)

temp = fft(subBinResponse,[],5);
subBinPh = angle(temp(:,:,:,:,2));
binMeanPh = squeeze(circ_mean(subBinPh));
binStdPh = squeeze(circ_std(subBinPh));
subBinPhDiff = squeeze(circ_dist(subBinPh(:,:,:,1),subBinPh(:,:,:,2)));
binMeanPhDiff = squeeze(mean(subBinPhDiff));
binStdPhDiff = squeeze(std(subBinPhDiff));

subBinDiff = squeeze(subBinAmp(:,:,:,1) - subBinAmp(:,:,:,2));%(iSub,iRoi,ibin)
binMeanDiff = squeeze(mean(subBinDiff));
binStdDiff = squeeze(std(subBinDiff));
binResponse = squeeze(mean(subBinResponse,1)); % mean timecourse across subjects. binResponse(iRoi,ibin,rwd,timepoint)
binResponseAmp = squeeze(std(binResponse,0,4)); %amplitude of mean timecourse. binResponseAmp(iRoi,ibin,rwd)



binVarMean = squeeze(mean(subBinVar,1)); % mean variability across subjects. binStdMean(iRoi,ibin,rwd,timepoint)
binVarStd = squeeze(std(subBinVar,0,1)); % std of variability across subjects. binStdStd(iRoi,ibin,rwd,timepoint)
binSubVarDiff = squeeze(subBinVar(:,:,:,2) - subBinVar(:,:,:,1));
binVarDiffMean = squeeze(mean(binSubVarDiff));
binVarDiffStd = squeeze(std(binSubVarDiff));

binStdVarMean = squeeze(mean(subBinStdVar));%mean trial std variability
binStdVarStd = squeeze(std(subBinStdVar));%std across subjects, of trial std variability
binSubStdVarDiff = squeeze(subBinStdVar(:,:,:,2) - subBinStdVar(:,:,:,1));
binStdVarDiffMean = squeeze(mean(binSubStdVarDiff));
binStdVarDiffStd = squeeze(std(binSubStdVarDiff));

% binStdMean = squeeze(mean(subBinStdMean));%mean trial std mean

binFftAmpVarMean = squeeze(mean(subBinFftAmpVar,1)); % mean variability across subjects. binStdMean(iRoi,ibin,rwd,timepoint)
binFftAmpVarStd = squeeze(std(subBinFftAmpVar,0,1)); % std of variability across subjects. binStdStd(iRoi,ibin,rwd,timepoint)
binSubFftAmpVarDiff = squeeze(subBinFftAmpVar(:,:,:,2) - subBinFftAmpVar(:,:,:,1));
binFftAmpVarDiffMean = squeeze(mean(binSubFftAmpVarDiff));
binFftAmpVarDiffStd = squeeze(std(binSubFftAmpVarDiff));

binFftPhVarMean = squeeze(mean(subBinFftPhVar,1)); % mean variability across subjects. binStdMean(iRoi,ibin,rwd,timepoint)
binFftPhVarStd = squeeze(std(subBinFftPhVar,0,1)); % std of variability across subjects. binStdStd(iRoi,ibin,rwd,timepoint)
binSubFftPhVarDiff = squeeze(subBinFftPhVar(:,:,:,2) - subBinFftPhVar(:,:,:,1));
binFftPhVarDiffMean = squeeze(mean(binSubFftPhVarDiff));
binFftPhVarDiffStd = squeeze(std(binSubFftPhVarDiff));




%%
i=1;
figure(i)
clf
rows=2;
cols = 9;
subplots = {1:3, 4:6, 7 , 9, cols+1:cols+3, cols+4:cols+6, cols+7 , cols+9};

for r= 1:length(ROIs)
    %amplitude
    subplot(rows,cols,subplots{1});
    iRoi=ROIs(r);
    for rwd=1:2
        dsErrorsurface(binCenters, squeeze(binMeanAmp(iRoi,:,rwd)), squeeze(binStdAmp(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(binCenters, squeeze(binMeanAmp(iRoi,:,rwd)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
    end
    subplot(rows,cols,subplots{2});
    dsErrorsurface(binCenters, binMeanDiff(iRoi,:), binStdDiff(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
    hold on
    plot(binCenters, squeeze(binMeanAmp(iRoi,:,1) -binMeanAmp(iRoi,:,2)),'k.-','linewidth',linewidth,'markersize',markersize);
    hline(0);
    
    %latency
%     binMeanPh = mod(binMeanPh,2*pi);
    binMeanPh = binMeanPh - pi/2;
    binMeanPh = mod((binMeanPh + pi),2*pi)-pi;
    subplot(rows,cols,cols+subplots{1});
    iRoi=ROIs(r);
    for rwd=1:2
        dsErrorsurface(binCenters, squeeze(binMeanPh(iRoi,:,rwd)), squeeze(binStdPh(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(binCenters, squeeze(binMeanPh(iRoi,:,rwd)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
    end
    subplot(rows,cols,cols+subplots{2});
    dsErrorsurface(binCenters, binMeanPhDiff(iRoi,:), binStdPhDiff(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
    hold on
    plot(binCenters, binMeanPhDiff(iRoi,:),'k.-','linewidth',linewidth,'markersize',markersize);
    hline(0);
    
end

%% bar plot of subjects' amplitudes
lineLength = 0.2;
lineWidth = 2;
markerSize = 10;
subplot(rows,cols,subplots{3})
for r = 1:length(ROIs)
    iROI = ROIs(r);
    clear smallSubAmp subPh
    for i=1:length(goodSubs)
        iSub = goodSubs(i);
        for rwd=1:2
            smallSubAmp(i,rwd) = std(squeeze(subResponse(iSub,iRoi,rwd,:)));
            f=angle(fft(squeeze(subResponse(iSub,iRoi,rwd,:))));
            subPh(i,rwd) = f(2);
        end
    end
    smallSubDiff = smallSubAmp(:,1) - smallSubAmp(:,2);
    minSubDiff = min(smallSubDiff);
    maxSubDiff= max(smallSubDiff);
    subjects = size(smallSubAmp,1);
    rewards = size(smallSubAmp,2);
    [rwdNum, subNum] = meshgrid(1:rewards, 1:subjects);
    scatterCmap=cool;
%     scatterCmap = scatterCmap(end:-1:1,:);%invert color map
    colormap(scatterCmap);

    l = size(scatterCmap,1);
    for i=1:length(goodSubs)
        subColor(i,:) = scatterCmap(1+floor((smallSubDiff(i) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
        plot(1:rewards,smallSubAmp(i,:),'Color',subColor(i,:),'linewidth',linewidth);
        hold on
    end
    scatter(rwdNum(:),smallSubAmp(:),markerSize, [subColor; subColor] ,'filled');   
    % MEAN
    for i=1:rewards
        scatter(i,mean(smallSubAmp(:,i)),markerSize*4,[0 0 0],'filled');
        plot(1:rewards,mean(smallSubAmp),'Color',[0 0 0],'linewidth', 2*linewidth);
    end
end

%% bar plot of subjects' latency
subplot(rows,cols,subplots{3+4})
for r = 1:length(ROIs)
    iROI = ROIs(r);
    subPhDiff = circ_dist(subPh(:,:),subPh(:,1));
    colormap(scatterCmap);
    for i=1:length(goodSubs)
        plot(1:rewards,subPhDiff(i,:),'Color',subColor(i,:),'linewidth',linewidth);
        hold on
    end
    scatter(rwdNum(:),subPhDiff(:),markerSize, [subColor; subColor] ,'filled');   
    % MEAN
    for i=1:rewards
        scatter(i,mean(subPhDiff(:,i)),markerSize*4,[0 0 0],'filled');
        plot(1:rewards,mean(subPhDiff),'Color',[0 0 0],'linewidth', 2*linewidth);
    end
%     ylim([-pi pi]);
end

%% Formatting
for isubplot=[1:3 5:7]%1:length(subplots)
    subplot(rows,cols,subplots{isubplot});
    switch isubplot
        case 1
            ylabel('response amplitude (std)');
        case 2
            ylabel('\Delta response amplitude (std)');
        case 3
             ylabel('response amplitude (std)');
        case 5
            ylabel('response timing (rad)');
        case 6
            ylabel('\Delta response timing (rad)');
        case 7
            ylabel('relative response latency (rad)');
    end
             
if isubplot==3
    xlabel('reward');
    drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -62/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
        'xTick',[1 2], 'xTickLabel', {'high','low'},...
        'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',7,...
        'yAxisMajorTickLen',-4/32);
elseif isubplot==7
    xlabel('reward');
    drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -52/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
        'xTick',[1 2], 'xTickLabel', {'high','low'},...
        'yTick',[-pi 0 pi], 'yTickLabel', {'-\pi','0','\pi'},...
        'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',7,...
        'yAxisMajorTickLen',-4/32);
    else
    
    xlabel('eccentricity (deg)');
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
        'labelFontSize',7);
    axis square
end
end
% set(gcf,'position',[10 10 25 	12]);
set(gcf,'position',[10 10 28 	14]);
% set(gcf,'position',[10 10 40 2*10]);
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig5_' ConcatProjStr '.pdf']);




%%
rows=3;
cols = 9;
subplots = {1:3, 4:6, 7 , 9};
iRoi = ROIs(1);

i=i+1;
figure(i)
clf

%% timepoint variability
subplot(rows,cols,subplots{1});
for rwd=1:2
    dsErrorsurface(binCenters, squeeze(binVarMean(iRoi,:,rwd)), squeeze(binVarStd(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot(binCenters, squeeze(binVarMean(iRoi,:,rwd)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
end
ylabel('mean variability (std)');
subplot(rows,cols,subplots{2});
dsErrorsurface(binCenters, binVarDiffMean(iRoi,:), binVarDiffStd(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
hold on
plot(binCenters, binVarDiffMean(iRoi,:),'k.-','linewidth',linewidth,'markersize',markersize);
ylabel('\Delta mean variability (std)');
hline(0);
%% STD amplitude variability
subplot(rows,cols,cols+subplots{1});
for rwd=1:2
    dsErrorsurface(binCenters, squeeze(binStdVarMean(iRoi,:,rwd)), squeeze(binStdVarStd(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot(binCenters, squeeze(binStdVarMean(iRoi,:,rwd)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
end
ylabel('amplitude variability (std)');
subplot(rows,cols,cols+subplots{2});
dsErrorsurface(binCenters, binStdVarDiffMean(iRoi,:), binStdVarDiffStd(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
hold on
plot(binCenters, binStdVarDiffMean(iRoi,:),'k.-','linewidth',linewidth,'markersize',markersize);
ylabel('\Delta amplitude variability (std)');
hline(0);
%% FFT amplitude variability
% subplot(rows,cols,cols+subplots{1});
% for rwd=1:2
%     dsErrorsurface(binCenters, squeeze(binFftAmpVarMean(iRoi,:,rwd)), squeeze(binFftAmpVarStd(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
%     hold on
%     plot(binCenters, squeeze(binFftAmpVarMean(iRoi,:,rwd)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
% end
% ylabel('amplitude variability (std)');
% subplot(rows,cols,cols+subplots{2});
% dsErrorsurface(binCenters, binFftAmpVarDiffMean(iRoi,:), binFftAmpVarDiffStd(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
% hold on
% plot(binCenters, binFftAmpVarDiffMean(iRoi,:),'k.-','linewidth',linewidth,'markersize',markersize);
% ylabel('\Delta amplitude variability (std)');
% hline(0);
%% FFT phase variability

subplot(rows,cols,2*cols+subplots{1});
for rwd=1:2
    dsErrorsurface(binCenters, squeeze(binFftPhVarMean(iRoi,:,rwd)), squeeze(binFftPhVarStd(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot(binCenters, squeeze(binFftPhVarMean(iRoi,:,rwd)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
end
ylabel('temporal variability (std)');
subplot(rows,cols,2*cols+subplots{2});
dsErrorsurface(binCenters, binFftPhVarDiffMean(iRoi,:), binFftPhVarDiffStd(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
hold on
plot(binCenters, binFftPhVarDiffMean(iRoi,:),'k.-','linewidth',linewidth,'markersize',markersize);
ylabel('\Delta temporal variability (std)');
hline(0);



%% get subject colors for bar plot
clear smallSubAmp
for i=1:length(goodSubs)
    iSub = goodSubs(i);
    for rwd=1:2
        smallSubAmp(i,rwd) = std(squeeze(subResponse(iSub,iRoi,rwd,:)));
    end
end
smallSubDiff = smallSubAmp(:,1) - smallSubAmp(:,2);
minSubDiff = min(smallSubDiff);
maxSubDiff= max(smallSubDiff);
subjects = size(smallSubAmp,1);
rewards = size(smallSubAmp,2);
[rwdNum, subNum] = meshgrid(1:rewards, 1:subjects);
scatterCmap=cool;
%     scatterCmap = scatterCmap(end:-1:1,:);%invert color map
colormap(scatterCmap);
l = size(scatterCmap,1);
for i=1:length(goodSubs)
    subColor(i,:) = scatterCmap(1+floor((smallSubDiff(i) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
end

%% compute subject variability for whole ROI
clear subVar subMeanVar roiBinFftAmp roiBinFftPh roiBinFftOther subBinFftAmpVar subBinFftPhVar subBinFftOtherVar
for iSub = 1:length(goodSubs)%length(subdirs)
    for rwd=1:2
        reshapedTrials = reshape(subTrialResponse{goodSubs(iSub),iRoi,rwd},10,[]);
        %trial-by-trial variability per subject
        subVar(iSub,rwd,:) = std(reshapedTrials,0,2);
        subMeanVar(iSub,rwd) = mean(squeeze(subVar(iSub,rwd,:)));
        temp = fft(reshapedTrials);
        roiBinFftAmp = abs(temp(2,:));
        roiBinFftPh = angle(temp(2,:));
        roiBinFftOther = sum(abs(temp([1 3:end],:)));
        subBinFftAmpVar(iSub,rwd) = std(roiBinFftAmp);
        subBinFftPhVar(iSub,rwd) = circ_std(roiBinFftPh');
        subBinFftOtherVar(iSub,rwd) = std(roiBinFftOther);
        
        singleTrialStd = std(reshapedTrials);
        subStdVar(iSub,rwd) = std(singleTrialStd);
        subStdMean(iSub,rwd) = mean(singleTrialStd);
    end
end

%% bar plot of subjects' timepoint variability
lineLength = 0.2;
lineWidth = 2;
markerSize = 10;

subplot(rows,cols,subplots{3})
for i=1:length(goodSubs)
    plot(1:rewards,subMeanVar(i,:),'Color',subColor(i,:),'linewidth',linewidth);
    hold on
end
scatter(rwdNum(:),subMeanVar(:),markerSize, [subColor; subColor] ,'filled');
% MEAN
for i=1:rewards
    scatter(i,mean(subMeanVar(:,i)),markerSize*4,[0 0 0],'filled');
    plot(1:rewards,mean(subMeanVar),'Color',[0 0 0],'linewidth', 2*linewidth);
end
ylabel('mean variability (std)');

%% bar plot of subjects' STD Amplitude variability
subplot(rows,cols,cols+subplots{3})
for i=1:length(goodSubs)
    plot(1:rewards,subStdVar(i,:),'Color',subColor(i,:),'linewidth',linewidth);
    hold on
end
scatter(rwdNum(:),subStdVar(:),markerSize, [subColor; subColor] ,'filled');
% MEAN
for i=1:rewards
    scatter(i,mean(subStdVar(:,i)),markerSize*4,[0 0 0],'filled');
    plot(1:rewards,mean(subStdVar),'Color',[0 0 0],'linewidth', 2*linewidth);
end
ylabel('amplitude variability (std)');

%% bar plot of subjects' FFT Amplitude variability
% subplot(rows,cols,cols+subplots{3})
% for i=1:length(goodSubs)
%     plot(1:rewards,subBinFftAmpVar(i,:),'Color',subColor(i,:),'linewidth',linewidth);
%     hold on
% end
% scatter(rwdNum(:),subBinFftAmpVar(:),markerSize, [subColor; subColor] ,'filled');
% % MEAN
% for i=1:rewards
%     scatter(i,mean(subBinFftAmpVar(:,i)),markerSize*4,[0 0 0],'filled','d');
%     plot(1:rewards,mean(subBinFftAmpVar),'Color',[0 0 0],'linewidth', 2*linewidth);
% end
% ylabel('amplitude variability (std)');

%% bar plot of subjects' FFT Phase variability
subplot(rows,cols,2*cols+subplots{3})
for i=1:length(goodSubs)
    plot(1:rewards,subBinFftPhVar(i,:),'Color',subColor(i,:),'linewidth',linewidth);
    hold on
end
scatter(rwdNum(:),subBinFftPhVar(:),markerSize, [subColor; subColor] ,'filled');
% MEAN
for i=1:rewards
    scatter(i,mean(subBinFftPhVar(:,i)),markerSize*4,[0 0 0],'filled');
    plot(1:rewards,mean(subBinFftPhVar),'Color',[0 0 0],'linewidth', 2*linewidth);
end
ylabel('temporal variability (std)');




%% formatting
for r=1:rows
    for isubplot=1:length(subplots)-1
        subplot(rows,cols,(r-1)*cols+subplots{isubplot})
%         ylabel('mean variability (std)');
        if isubplot<3
            xlabel('eccentricity');
%             if isubplot==2 
%                 ylabel('\Delta mean variability (std)');
%             else
%                 ylabel('mean variability (std)');
%             end
            
            drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
                'labelFontSize',7);
            axis square
        else
            xlabel('reward');
            drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -50/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
                'xTick',[1 2], 'xTickLabel', {'high','low'},...
                'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',7,...
                'yAxisMajorTickLen',-4/32);
        end
    end
end


set(gcf,'position',[10 10 25 	24]);

print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig5_variability_' ConcatProjStr '.pdf']);


%%

[rmat pmat] = corr([smallSubDiff,subMeanVar(:,2)-subMeanVar(:,1),subBinFftAmpVar(:,2)-subBinFftAmpVar(:,1),...
    subBinFftPhVar(:,2)-subBinFftPhVar(:,1), subStdMean(:,2)-subStdMean(:,1),subStdVar(:,2)-subStdVar(:,1)]);
lbls = {'std amp','timepoint var','FFT amp var','FFT ph var','std amp mean','std amp var'};

pmat


