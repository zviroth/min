close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

% randSeed = rng;
load(fullfile(dataFolder,'randSeed.mat'),'randSeed');
rng(randSeed);
% save(fullfile(dataFolder,'randSeed.mat'),'randSeed');


onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    
    onlyCorrectString = '_incorrect';
elseif onlyCorrect==0
    onlyCorrectString = '_validresponse';
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
    'subRoiRuns','runMeanFFT',...
    'allVoxTrialResponse','allVoxTaskPhase','allVoxTaskAmp','allVoxTaskCo',...
    'voxTrials','voxGoodTrials','meanVoxTrial',...
    'maxRT');

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
fontsize=9;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);

goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% goodSubs = [1 3 5:10 13]; %excluding 4 subjects with highest amplitude difference

% ROIs = [1:4];
% ROIs = [1 2];
rwdString = {'H','L'};
%% ANALYZE GLOBAL MEAN SIGNAL
for iSub = 1:length(goodSubs)
    for rwd=1:2
         globalTrials = reshape(globalMean{iSub,rwd},trialLength,[]);
         meanGlobalTrial(iSub,rwd,:) = mean(globalTrials');
         varGlobalTrial(iSub,rwd,:) = std(globalTrials');
    end
end

for rwd=1:2
    plot(squeeze(mean(meanGlobalTrial(:,rwd,:))),'color',plotColors{rwd});
    hold on
end


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

for iSub = 1:length(goodSubs)%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(roiNames)
        for rwd=1:2
%             reshapedTrials = reshape(subTrialResponse{goodSubs(iSub),iRoi,rwd},10,[]);
            reshapedTrials = subTrialResponse{goodSubs(iSub),iRoi,rwd};
           temp = fft(reshapedTrials);
           roiFftAmp = abs(temp(2,:));
           roiFftBaseline = abs(temp(1,:));

           subTrialResponse{goodSubs(iSub),iRoi,rwd} = reshapedTrials;
    
        %trial-by-trial variability per subject
           subTimepointStd(iSub,iRoi,rwd,:) = std(reshapedTrials,0,2);
           temp = fft(reshapedTrials);
           roiFftAmp = abs(temp(2,:));
           roiFftPh = angle(temp(2,:));
           trialRoiFftPhVec{iSub,iRoi,rwd} = roiFftPh;
           trialRoiFftAmpVec{iSub,iRoi,rwd} = roiFftAmp;
           trialRoiAmpVec{iSub,iRoi,rwd} = std(reshapedTrials);
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

%         if rwd==1
%             firstTrial(rwd)=1;
%         else %rwd==2
%             firstTrial(rwd)=numTrials(iSub,1)+1;
%         end
    end
    firstTrial(1)=1;
    firstTrial(2)=numTrials(iSub,1)+1;
    
    for iRoi=ROIs%length(roiNames)
        roiTrials{iRoi} = [subTrialResponse{goodSubs(iSub),iRoi,1} subTrialResponse{goodSubs(iSub),iRoi,2}];
    end
    subRT = [trialRTvec{iSub,1}; trialRTvec{iSub,2}];
    
   
    for iRoi= ROIs%length(roiNames)
        %difference in amplitude
        for rwd = 1:2
            realSubMeanTC(iSub,iRoi,rwd,:) = mean(subTrialResponse{goodSubs(iSub),iRoi,rwd},2);
            temp = fft(squeeze(realSubMeanTC(iSub,iRoi,rwd,:)));
            realSubMeanFftAmp(iSub,iRoi,rwd) = abs(temp(2));
            realSubMeanFftPh(iSub,iRoi,rwd) = angle(temp(2));
        end
        realSubDiff(iSub,iRoi) = std(realSubMeanTC(iSub,iRoi,1,:)) - std(realSubMeanTC(iSub,iRoi,2,:));%difference between high std and low std

        %difference in latency
        realSubFftAmpDiff(iSub,iRoi) = realSubMeanFftAmp(iSub,iRoi,1) - realSubMeanFftAmp(iSub,iRoi,2);
        realSubFftPhDiff(iSub,iRoi) = circ_dist(realSubMeanFftPh(iSub,iRoi,2), realSubMeanFftPh(iSub,iRoi,1));
        
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

%RT variability
realRTvar = mean(realSubRTvar);
realRTvarDiff = realRTvar(2) - realRTvar(1);
%RT
realMeanRT = mean(realSubRT);
realStdRT = std(realSubRT);
realRTdiff = realMeanRT(2) - realMeanRT(1);

%%
% histbins=1000;
% rows=length(roiNames);
% cols = length(goodSubs);


%%

realSubRTvarDiff = realSubRTvar(:,2) - realSubRTvar(:,1);
realSubRTdiff = realSubRT(:,2)-realSubRT(:,1);

%% correlate trial-by-trial RT with response amplitude/timing
%correlate RT with BOLD timing


for iRoi= ROIs
    for rwd=1:2
        allRTvec{iRoi,rwd}=[];
        allPhVec{iRoi,rwd}=[];
        allStdVec{iRoi,rwd}=[];
        allAmpVec{iRoi,rwd}=[];
        for iSub=1:length(goodSubs)
            [rtBoldPhCorr(iSub,iRoi,rwd) rtBoldPhCorrPval(iSub,iRoi,rwd)]=  corr(trialRoiFftPhVec{iSub,iRoi,rwd}', trialRTvec{iSub,rwd});
            [rtBoldPhShiftCorr(iSub,iRoi,rwd) rtBoldPhShiftCorrPval(iSub,iRoi,rwd)]=  corr(circshift(trialRoiFftPhVec{iSub,iRoi,rwd}',-1), trialRTvec{iSub,rwd});
            
            [rtBoldAmpCorr(iSub,iRoi,rwd) rtBoldAmpCorrPval(iSub,iRoi,rwd)]=  corr(trialRoiFftAmpVec{iSub,iRoi,rwd}', trialRTvec{iSub,rwd});
            [rtBoldAmpShiftCorr(iSub,iRoi,rwd) rtBoldAmpShiftCorrPval(iSub,iRoi,rwd)]=  corr(circshift(trialRoiFftAmpVec{iSub,iRoi,rwd}',-1), trialRTvec{iSub,rwd});
            
            [rtBoldStdCorr(iSub,iRoi,rwd) rtBoldStdCorrPval(iSub,iRoi,rwd)]=  corr(trialRoiAmpVec{iSub,iRoi,rwd}', trialRTvec{iSub,rwd});
            [rtBoldStdShiftCorr(iSub,iRoi,rwd) rtBoldStdShiftCorrPval(iSub,iRoi,rwd)]=  corr(circshift(trialRoiAmpVec{iSub,iRoi,rwd}',-1), trialRTvec{iSub,rwd});
            
            
            subTrialMat = [trialRoiFftPhVec{iSub,iRoi,rwd}', trialRoiFftAmpVec{iSub,iRoi,rwd}',trialRoiAmpVec{iSub,iRoi,rwd}',-trialRTvec{iSub,rwd}];
            [subTrialCorr(iSub,iRoi,rwd,:,:) subTrialCorrPval(iSub,iRoi,rwd,:,:)] = corr(subTrialMat);
            
            %pool across subjects
%             allRTvec{iRoi,rwd} = [allRTvec{iRoi,rwd}; zscore(trialRTvec{iSub,rwd})]; 
%             allPhVec{iRoi,rwd}= [allPhVec{iRoi,rwd}; zscore(trialRoiFftPhVec{iSub,iRoi,rwd})'];
%             allStdVec{iRoi,rwd}=[allStdVec{iRoi,rwd}; zscore(trialRoiAmpVec{iSub,iRoi,rwd})'];
%             allAmpVec{iRoi,rwd}=[allAmpVec{iRoi,rwd}; zscore(trialRoiFftAmpVec{iSub,iRoi,rwd})'];
            
            allRTvec{iRoi,rwd} = [allRTvec{iRoi,rwd}; trialRTvec{iSub,rwd}]; 
            allPhVec{iRoi,rwd}= [allPhVec{iRoi,rwd}; trialRoiFftPhVec{iSub,iRoi,rwd}'];
            allStdVec{iRoi,rwd}=[allStdVec{iRoi,rwd}; trialRoiAmpVec{iSub,iRoi,rwd}'];
            allAmpVec{iRoi,rwd}=[allAmpVec{iRoi,rwd}; trialRoiFftAmpVec{iSub,iRoi,rwd}'];
        end
        rwdTrialMat{iRoi, rwd} = [allPhVec{iRoi,rwd}'; allAmpVec{iRoi,rwd}'; allStdVec{iRoi,rwd}'; -allRTvec{iRoi,rwd}'];
        [roiTrialCorr(iRoi,rwd,:,:) roiTrialCorrPval(iRoi,rwd,:,:)] = corr(rwdTrialMat{iRoi, rwd}');
    end
    trialMat{iRoi} = [rwdTrialMat{iRoi, 1} rwdTrialMat{iRoi, 2}]; 
    [roiTrialCorr(iRoi,3,:,:) roiTrialCorrPval(iRoi,3,:,:)] = corr(trialMat{iRoi}');
    [circLinCorr(iRoi) circLinCorrPval(iRoi)] = circ_corrcc([allPhVec{iRoi,1}' allPhVec{iRoi,2}'], -[allRTvec{iRoi,1}' allRTvec{iRoi,2}']);
end
trialLabels = {'FFT ph','FFT amp','STD amp','-RT'};


%permutation test of pooled correlations
nfactors = size(rwdTrialMat{iRoi,rwd},1);
ntotaltrials = size(rwdTrialMat{iRoi,rwd},2);



%% mean over single subject correlations
for iRoi=ROIs
    for rwd=1:2
        meanSubTrialCorr(iRoi,rwd,:,:) = squeeze(mean(subTrialCorr(:,iRoi,rwd,:,:)));
    end
end

%%
mean(subPropCorrect)
std(subPropCorrect)
mean(realSubRT)
std(realSubRT)

%%
i=0;

i=i+1; figure(i); clf
iRoi=2;
rows=2;
cols = 9;
% subplots = {1:3, 4:6, 7 , 9, cols+1:cols+3, cols+4:cols+6, cols+7 , cols+9};
subplots = {1 3 5 7:9};
lineLength = 0.2;
lineWidth = 2;
markerSize = 10;


for iSub=1:length(goodSubs)
    for rwd=1:2
        smallSubAmp(iSub,rwd) = std(squeeze(subResponse(goodSubs(iSub),iRoi,rwd,:)));
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

% RT
l = size(scatterCmap,1);
% percent correct
subplot(rows,cols,subplots{1})
for iSub=1:length(goodSubs)
    subColor(iSub,:) = scatterCmap(1+floor((smallSubDiff(iSub) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
    plot(1:rewards,subPropCorrect(iSub,:),'Color',subColor(iSub,:),'linewidth',linewidth);
    hold on
end
scatter(rwdNum(:),subPropCorrect(:),markerSize, [subColor; subColor] ,'filled');
% MEAN
for irwd=1:rewards
    scatter(irwd,mean(subPropCorrect(:,irwd)),markerSize*4,[0 0 0],'filled');
    plot(1:rewards,mean(subPropCorrect),'Color',[0 0 0],'linewidth', 2*linewidth);
end

%RT
subplot(rows,cols,subplots{2})
for i=1:length(goodSubs)
    plot(1:rewards,realSubRT(i,:),'Color',subColor(i,:),'linewidth',linewidth);
    hold on
end
scatter(rwdNum(:),realSubRT(:),markerSize, [subColor; subColor] ,'filled');
% MEAN
for i=1:rewards
    scatter(i,mean(realSubRT(:,i)),markerSize*4,[0 0 0],'filled');
    plot(1:rewards,mean(realSubRT),'Color',[0 0 0],'linewidth', 2*linewidth);
end

% RT variability
subplot(rows,cols,subplots{3})
for iSub=1:length(goodSubs)
%     subColor(i,:) = scatterCmap(1+floor((smallSubDiff(i) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
    plot(1:rewards,realSubRTvar(iSub,:),'Color',subColor(iSub,:),'linewidth',linewidth);
    hold on
end
scatter(rwdNum(:),realSubRTvar(:),markerSize, [subColor; subColor] ,'filled');
% MEAN
for irwd=1:rewards
    scatter(irwd,mean(realSubRTvar(:,irwd)),markerSize*4,[0 0 0],'filled');
    plot(1:rewards,mean(realSubRTvar),'Color',[0 0 0],'linewidth', 2*linewidth);
end

% RT vs. latency
subplot(rows,cols,subplots{4})
scatter(allPhVec{iRoi,rwd},allRTvec{iRoi,rwd},markerSize/2,[0 0 0]);
xlabel('fMRI response latency (rad)');
ylabel('RT (ms)');
drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -16/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
    'labelFontSize',fontsize,'xAxisMin', -pi ,'xAxisMax', pi,'xTick',[-pi pi],'xTickLabel', {'-\pi','\pi'});
axis square
%
for isubplot=1:3
    subplot(rows,cols,subplots{isubplot})
    switch isubplot
        case 1
            ylabel('correct (%)');
        case 2
            ylabel('RT (ms)');
        case 3
            ylabel('RT variability (ms)');
    end 
    xlabel('reward');
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -80/64, 'xAxisMargin', 4/64, 'yAxisMargin', 4/64,'xAxisMinMaxSetByTicks',1,...
        'xTick',[1 2], 'xTickLabel', {'high','low'},...
        'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',fontsize,...
        'yAxisMajorTickLen',-4/32);
end

colorSubFolders = subFolders;
colorGoodSubs = goodSubs;
% save(fullfile(dataFolder,'subColors.mat'),'subColor', 'scatterCmap', 'colorSubFolders', 'colorGoodSubs');

set(gcf,'position',[10 10 21 	14]);
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig_RT.pdf']);


%%
for iSub=1:length(goodSubs)
    for rwd=1:2
%         trialResponse{goodSubs(iSub),rwd}(iRun,itrial);
        percentAbort(iSub,rwd) = 1 - sum(trialResponse{goodSubs(iSub),rwd}(:)>0) ./ length(trialResponse{goodSubs(iSub),rwd}(:));
    end
end