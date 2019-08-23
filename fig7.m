close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

nperms=10000;
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
load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'globalMean','regressBetasGlobal','runRwd',...
    'subRoiRuns','runMeanFFT');
% trialLength=10;
clear subMeanResponse trialStd meanTrialStd runStd meanRunStd subRwdStd trialFFTamp trialFFTphase meanTrialFFTamp meanTrialFFTphase 
clear stdTrialFFTphase runFFTamp runFFTphase meanRunFFTamp meanRunFFTphase stdRunFFTphase subFFTamp subFFTphase
clear trialMinMax meanTrialMinMax runMinMax meanRunMinMax subRwdMinMax legendCellArray
clear stdRunMax stdTrialMax
clear groupLabels pVal
clear yError yMin yMax permDiff meanPermDiff
clear permSubDiff realSubDiff permSubResp realSubResp meanPermTC permSubMeanTC realSubMeanTC meanRealTC
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
ROIs = [1:4];
% ROIs = [1 2];


i=0;
%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;



%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd

% goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% goodSubs = [1 3 5:10 13]; %excluding 4 subjects with highest amplitude difference


for iSub = 1:length(goodSubs)%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(roiNames)
        for rwd=1:2
%             reshapedTrials = reshape(subTrialResponse{goodSubs(iSub),iRoi,rwd},10,[]);
            reshapedTrials = subTrialResponse{goodSubs(iSub),iRoi,rwd};
           temp = fft(reshapedTrials);
           roiFftAmp = abs(temp(2,:));
           roiFftBaseline = abs(temp(1,:));
           
           %remove FFT component 1
%            reshapedTrials = reshapedTrials - repmat(mean(reshapedTrials),10,1);

           %remove FFT component 2
%            reshapedTrials = reshapedTrials./repmat(roiBinFftAmp,10,1);
%            reshapedTrials = (reshapedTrials - repmat(mean(reshapedTrials),10,1))./repmat(roiFftAmp,10,1) + repmat(mean(reshapedTrials),10,1);
           
           subTrialResponse{goodSubs(iSub),iRoi,rwd} = reshapedTrials;
    
        %trial-by-trial variability per subject
           subTimepointStd(iSub,iRoi,rwd,:) = std(reshapedTrials,0,2);
           temp = fft(reshapedTrials);
           roiFftAmp = abs(temp(2,:));
           roiFftPh = angle(temp(2,:));
           roiFftOther = sum(abs(temp([1 3:end],:)));
%            roiBinFftBaseline = abs(temp(1,:));
           roiFftBaseline = mean(reshapedTrials);
           subFftAmpStd(iSub,iRoi,rwd) = std(roiFftAmp);
           subFftMeanAmp(iSub,iRoi,rwd) = mean(roiFftAmp);
           subFftMeanPh(iSub,iRoi,rwd) = circ_mean(roiFftPh');
           subMeanStd(iSub,iRoi,rwd) = mean(std(reshapedTrials));
           subStdVar(iSub,iRoi,rwd) = std(std(reshapedTrials));%variability of amplitude measured by std
           subMedianStd(iSub,iRoi,rwd) = median(std(reshapedTrials));
           
           subFftPhStd(iSub,iRoi,rwd) = circ_std(roiFftPh');
           subFftOtherStd(iSub,iRoi,rwd) = std(roiFftOther);
           subFftBaselineStd(iSub,iRoi,rwd) = std(roiFftBaseline);

            if iSub==goodSubs(1)
                allTrials{iRoi,rwd} = reshapedTrials;
            else
                allTrials{iRoi,rwd} = [allTrials{iRoi,rwd} reshapedTrials];
            end
        end
     end
end


for rwd=1:2
    for iRoi = 1:length(roiNames)
        meanResponse(iRoi,rwd,:) = mean(subResponse(goodSubs,iRoi,rwd,:));
        stdResponse(iRoi,rwd,:) = std(subResponse(goodSubs,iRoi,rwd,:));
%         plot(squeeze(meanResponse(rwd,iRoi,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
%         hold on
        meanTimepointStd(iRoi, rwd,:) = mean(subTimepointStd(:,iRoi,rwd,:));%mean trial-by-trial variability
        stdTimepointStd(iRoi, rwd,:) = std(subTimepointStd(:,iRoi,rwd,:));%std of trial-by-trial variability
    end
end


%% PERMUTATIONS
tic

for iSub = 1:length(goodSubs)%length(subdirs)
    for rwd=1:2

        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!
        temp = trialResponse{goodSubs(iSub),rwd}(:,2:end-1);
%         temp = trialResponse{goodSubs(iSub),rwd}(:);
        trialResponseVec = temp(:);
%         temp = trialRT{goodSubs(iSub),rwd}(:);
        temp = trialRT{goodSubs(iSub),rwd}(:,2:end-1);
        trialRTvec{iSub,rwd} = temp(trialResponseVec>0);%only trials with a response, we're assuming onlyCorrect==0
        temp = trialCorrectness{goodSubs(iSub),rwd}(:,1:end);%include the first trial, and the last as well
        trialCorrectVec{iSub,rwd} = temp(:); %include trials with no response
        subPropCorrect(iSub,rwd) = sum(trialCorrectVec{iSub,rwd})/length(trialCorrectVec{iSub,rwd});

    end
    firstTrial(1)=1;
    firstTrial(2)=numTrials(iSub,1)+1;
    
    for iRoi=ROIs%length(roiNames)
        roiTrials{iRoi} = [subTrialResponse{goodSubs(iSub),iRoi,1} subTrialResponse{goodSubs(iSub),iRoi,2}];
    end
    subRT = [trialRTvec{iSub,1}; trialRTvec{iSub,2}];

    for iRoi= ROIs%length(roiNames)
        for rwd = 1:2
            realSubMeanTC(iSub,iRoi,rwd,:) = mean(subTrialResponse{goodSubs(iSub),iRoi,rwd},2);
        end
        realSubDiff(iSub,iRoi) = std(realSubMeanTC(iSub,iRoi,1,:)) - std(realSubMeanTC(iSub,iRoi,2,:));%difference between high std and low std

        %FFT variability
        realSubAmpStdDiff(iSub,iRoi) = subFftAmpStd(iSub,iRoi,2) - subFftAmpStd(iSub,iRoi,1);
        realSubPhStdDiff(iSub,iRoi) = subFftPhStd(iSub,iRoi,2) - subFftPhStd(iSub,iRoi,1);
        realSubOtherStdDiff(iSub,iRoi) = subFftOtherStd(iSub,iRoi,2) - subFftOtherStd(iSub,iRoi,1);
        realSubBaselineStdDiff(iSub,iRoi) = subFftBaselineStd(iSub,iRoi,2) - subFftBaselineStd(iSub,iRoi,1);
        
        realSubMeanFftAmpDiff(iSub,iRoi) = subFftMeanAmp(iSub,iRoi,2) - subFftMeanAmp(iSub,iRoi,1);
        realSubMeanFftPhDiff(iSub,iRoi) = squeeze(circ_dist(subFftMeanPh(iSub,iRoi,2),subFftMeanPh(iSub,iRoi,1)));
        realSubMeanStdDiff(iSub,iRoi) = subMeanStd(iSub,iRoi,2) - subMeanStd(iSub,iRoi,1);
        realSubMedianStdDiff(iSub,iRoi) = subMedianStd(iSub,iRoi,2) - subMedianStd(iSub,iRoi,1);
        
        realSubStdVarDiff(iSub,iRoi) = subStdVar(iSub,iRoi,2) - subStdVar(iSub,iRoi,1);

    end
    for rwd=1:2
        %RT variability
        realSubRTvar(iSub,rwd) = std(trialRTvec{iSub,rwd});
        %RT
        realSubRT(iSub,rwd) = mean(trialRTvec{iSub,rwd});
    end
    
end


meanRealAmpStdDiff = squeeze(mean(realSubAmpStdDiff));%mean over subjects
meanRealPhStdDiff = squeeze(mean(realSubPhStdDiff));%mean over subjects
meanRealOtherStdDiff = squeeze(mean(realSubOtherStdDiff));%mean over subjects
meanRealBaselineStdDiff = squeeze(mean(realSubBaselineStdDiff));%mean over subjects
meanRealFftAmpDiff = squeeze(mean(realSubMeanFftAmpDiff));%mean over subjects
meanRealFftPhDiff = squeeze(mean(realSubMeanFftPhDiff));%mean over subjects
meanRealStdVarDiff = squeeze(mean(realSubStdVarDiff));%mean over subjects

meanRealMeanStdDiff = squeeze(mean(realSubMeanStdDiff));%mean trials amplitude, mean over subjects
meanRealMedianStdDiff = squeeze(mean(realSubMedianStdDiff));%median trials amplitude, mean over subjects


for iRoi= ROIs
    %RFX
    meanRealDiff(iRoi) = mean(realSubDiff(:,iRoi));%average over subjects

    %FFX
    for rwd = 1:2
        meanRealTC(iRoi,rwd,:) = squeeze(mean(realSubMeanTC(:,iRoi,rwd,:)));%average timecourse over subjects
    end
    meanRealDiffStd(iRoi) = std(squeeze(meanRealTC(iRoi,1,:))) - std(squeeze(meanRealTC(iRoi,2,:)));


end

%RT variability

realRTvar = mean(realSubRTvar);
realRTvarDiff = realRTvar(2) - realRTvar(1);
%RT
realMeanRT = mean(realSubRT);
realStdRT = std(realSubRT);
realRTdiff = realMeanRT(2) - realMeanRT(1);

clear smallSubAmp
for iRoi = 1:length(roiNames)
    for s=1:length(goodSubs)
        iSub = goodSubs(s);
        for rwd=1:2
            smallSubAmp(s,iRoi,rwd) = std(squeeze(subResponse(iSub,iRoi,rwd,:)));
        end
    end
end


iRoi=2;

smallSubDiff = smallSubAmp(:,:,1) - smallSubAmp(:,:,2);
subStdDiff = mean(subTimepointStd(:,:,2,:),4) - mean(subTimepointStd(:,:,1,:),4);


realSubRTvarDiff = realSubRTvar(:,2) - realSubRTvar(:,1);

realSubRTdiff = realSubRT(:,2)-realSubRT(:,1);


%% CORRELATION MATRIX
iRoi=2;
[rmat pmat] = corr([smallSubDiff(:,iRoi),realSubMeanFftPhDiff(:,iRoi),subStdDiff(:,iRoi),realSubPhStdDiff(:,iRoi),realSubAmpStdDiff(:,iRoi),realSubStdVarDiff(:,iRoi),realSubRTvarDiff,realSubRTdiff])

 i=i+1; figure; clf
 rows=1; 
 cols=2;
subplot(rows,cols,1)
imagesc(rmat)
lbls = {'std amp','mean phase', 'timepoint var','FFT ph var','FFT amp var','std amp var','RT var','mean RT' };
yticks(1:length(lbls));
yticklabels(lbls);

 subplot(rows,cols,2)
 imagesc(pmat<0.05);
yticks(1:length(lbls));
yticklabels(lbls);


toc

%%
iRoi=2;
stdOfMeanTrial = mean(std(realSubMeanTC(:,iRoi,:,:),0,4));
meanOfTrialStd = mean(subMeanStd(:,iRoi,:));
meanOfTrialsFFTamp = mean(subFftMeanAmp(:,iRoi,:));
medianOfTrialsStd = subMedianStd(:,iRoi,:);

 temp = std(realSubMeanTC(:,iRoi,:,:),0,4);
 stdMeanTrialDiff = squeeze(temp(:,:,1) - temp(:,:,2));

 medianTrialStdDiff = squeeze(subMedianStd(:,iRoi,1) - subMedianStd(:,iRoi,2));
 
 %% MEAN SPECTRUM
iRoi=2;
 i=i+1; figure; clf
 TR=1.5;
 dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;
 for rwd=1:2
     runFFT = squeeze(mean(runMeanFFT(:,iRoi,rwd,:)));
     runFFTstd = squeeze(std(runMeanFFT(:,iRoi,rwd,:)));
     nframes = length(runFFT);
        %from plotMeanFourierAmp.m
        runFFT = runFFT / (nframes/2);
        frequencies = [0:nframes-1]/(nframes*TR);
        runFFT = runFFT(1:floor(nframes/2));
        frequencies = frequencies(1:floor(nframes/2));
        runFFTstd = runFFTstd / (nframes/2);
        runFFTstd = runFFTstd(1:floor(nframes/2));

        dsErrorsurface(frequencies, runFFT, runFFTstd./sqrt(size(runMeanFFT,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(frequencies, runFFT,'.-', 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);

 end
 
 %%
 %% Timepoint variability as function of time
subTimepointStdDiff = squeeze(subTimepointStd(:,:,1,:) - subTimepointStd(:,:,2,:));
 i=i+1; figure; clf
 cols=5;
 rows=ceil(length(goodSubs)/cols);
 for iSub=1:length(goodSubs)
     subplot(rows,cols,iSub)
     for rwd=1:2
         plot(squeeze(subTimepointStd(iSub,iRoi,rwd,:)),'color', plotColors{rwd},'linewidth',2);
         hold on
     end
     plot(squeeze(subTimepointStdDiff(iSub,iRoi,:)),'k','linewidth',2);
 end
 groupMeanVar = mean(subTimepointStd(:,iRoi,:,:));
 subplot(rows,cols,rows*cols)
 for rwd=1:2
     plot(squeeze(groupMeanVar(:,:,rwd,:)),'color', plotColors{rwd},'linewidth',3);
     hold on
 end
 plot(squeeze(groupMeanVar(:,:,2,:)-groupMeanVar(:,:,1,:)),'k','linewidth',2);
  title('mean timepoint variability');

%%
iRoi=2;
dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;


groupMeanVar = squeeze(mean(subTimepointStd(:,iRoi,:,:)));
groupStdVar = squeeze(std(subTimepointStd(:,iRoi,:,:)));

subVarDiff = squeeze(subTimepointStd(:,iRoi,2,:) - subTimepointStd(:,iRoi,1,:));
meanVarDiff = mean(subVarDiff);
stdVarDiff = std(subVarDiff);
i=i+1; figure; clf
subplots = {1:3, 4:6, 7 , 9};
rows=1;
cols = 9;
subplot(rows,cols,subplots{1})
for rwd=1:2
    dsErrorsurface(TR*(0:trialLength-1), groupMeanVar(rwd,:), groupStdVar(rwd,:)./sqrt(size(subTimepointStd,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot(TR*(0:trialLength-1),groupMeanVar(rwd,:), 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);
end
subplot(rows,cols,subplots{2})
dsErrorsurface(TR*(0:trialLength-1), meanVarDiff, stdVarDiff./sqrt(size(subTimepointStd,1)), [0 0 0],dsSurfaceAlpha);
hold on
plot(TR*(0:trialLength-1),meanVarDiff,'k-','linewidth',linewidth,'markersize',markersize);
hline(0);

for isubplot=1:2
    subplot(rows,cols,subplots{isubplot});
    xlabel('time (sec)');
    if isubplot==1
        ylabel('mean variability (std)');
    else
        ylabel('\Delta mean variability (std)');
    end
    drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',0,...
        'labelFontSize',7);
    axis square
end
set(gcf,'position',[10 10 25 	6]);
% set(gcf,'position',[10 10 40 2*10]);
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig7_' ConcatProjStr '.pdf']);
%%
mean(subPropCorrect)
std(subPropCorrect)
mean(realSubRT)
std(realSubRT)