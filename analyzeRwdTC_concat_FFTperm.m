close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj=1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

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
    'globalMean','regressBetasGlobal','runRwd');

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
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
% ROIs = [1 2 9 10];
% ROIs = [1 2];

%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;

%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd subStd permSubFftDiff

% goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% subjects. 0.0457  0.0140,     0.0553    0.0214


for rwd=1:2
    for iRoi = 1:length(roiNames)
        meanResponse(rwd,iRoi,:) = mean(subResponse(goodSubs,iRoi,rwd,:));
        stdResponse(rwd,iRoi,:) = std(subResponse(goodSubs,iRoi,rwd,:));
%         plot(squeeze(meanResponse(rwd,iRoi,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
%         hold on
    end
end


%% PERMUTATIONS
tic

for iSub = 1:length(goodSubs)%length(subdirs)
    for rwd=1:2

        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!
        numRuns(iSub,rwd) = size(subRoiRuns{goodSubs(iSub),1,rwd},2);
    end
    firstTrial(1)=1;
    firstTrial(2)=numTrials(iSub,1)+1;
    firstRun(1) = 1;
    firstRun(2)=numRuns(iSub,1)+1;
    
    for iRoi=ROIs%length(roiNames)
        roiTrials{iRoi} = [subTrialResponse{goodSubs(iSub),iRoi,1} subTrialResponse{goodSubs(iSub),iRoi,2}];
        roiRuns{iRoi} = [subRoiRuns{goodSubs(iSub),iRoi,1} subRoiRuns{goodSubs(iSub),iRoi,2}];
    end
    
    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));
        for iRoi= ROIs%length(roiNames)
            for rwd=1:2
                permTrials = roiTrials{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));% 10 timepoints X #trials
                permSubMeanTC(iSub,iRoi,rwd,p,:) = mean(permTrials,2);%average timecourse over trials
                permSubResp(iSub,iRoi,rwd,p) = std(squeeze(permSubMeanTC(iSub,iRoi,rwd,p,:)));%amplitude of task response
                
                permSubVar(iSub,iRoi,rwd,p,:) = std(permTrials,0,2);%variability timecourse over trials
            end
        end
        randOrder = randperm(numRuns(iSub,1)+numRuns(iSub,2));
        for iRoi= ROIs%length(roiNames)
            for rwd=1:2
                permRuns = roiRuns{iRoi}(:,randOrder(firstRun(rwd):firstRun(rwd)+numRuns(iSub,rwd)-1));
                permSubMeanRun(iSub,iRoi,rwd,p,:) = mean(permRuns,2);%average timecourse over runs
                temp = abs(fft(permSubMeanRun(iSub,iRoi,rwd,p,:)));
                permSubFFT(iSub,iRoi,rwd,p) = sum(temp(ntrials+1:ntrials:(1+end/2)));%amplitude of task response

            end
        end
    end
    
    for iRoi= ROIs%length(roiNames)
        for rwd = 1:2
            realSubMeanRun(iSub,iRoi,rwd,:) = mean(subRoiRuns{goodSubs(iSub),iRoi,rwd},2);
            realSubMeanTC(iSub,iRoi,rwd,:) = mean(subTrialResponse{goodSubs(iSub),iRoi,rwd},2);
            
           %trial-by-trial variability per subject
           subVar(iSub,iRoi,rwd,:) = std(subTrialResponse{goodSubs(iSub),iRoi,rwd},0,2);

        end
        realSubDiff(iSub,iRoi) = std(realSubMeanTC(iSub,iRoi,1,:)) - std(realSubMeanTC(iSub,iRoi,2,:));%difference between high std and low std
        permSubDiff(iSub,iRoi,:) = squeeze(permSubResp(iSub,iRoi,1,:) - permSubResp(iSub,iRoi,2,:));%difference between high std and low std
        pVal_sub(iSub, iRoi) = sum(permSubDiff(iSub,iRoi,:) > realSubDiff(iSub,iRoi))/nperms;
        
        temp = abs(fft(realSubMeanRun(iSub,iRoi,1,:))) - abs(fft(realSubMeanRun(iSub,iRoi,2,:)));
        realSubFftDiff(iSub,iRoi) = sum(temp(ntrials+1:ntrials:(1+end/2)));
        permSubFftDiff(iSub,iRoi,:) = permSubFFT(iSub,iRoi,1,:) - permSubFFT(iSub,iRoi,2,:);

    end
end
hypoth = pVal_sub<0.05;
% pVal_sub(:,ROIs);
for iRoi= ROIs
    %RFX
    meanPermDiff(iRoi,:) = mean(permSubDiff(:,iRoi,:));%average over subjects
    meanRealDiff(iRoi) = mean(realSubDiff(:,iRoi));%average over subjects
    pVal_rfx(iRoi) = sum(meanPermDiff(iRoi,:) > meanRealDiff(iRoi))/nperms;
    
    %FFT RFX
    meanPermFftDiff(iRoi,:) = mean(permSubFftDiff(:,iRoi,:));%average over subjects
    meanRealFftDiff(iRoi) = mean(realSubFftDiff(:,iRoi));%average over subjects
    pVal_fft_rfx(iRoi) = sum(meanPermFftDiff(iRoi,:) > meanRealFftDiff(iRoi))/nperms;
    
    %FFX
    for rwd = 1:2
        meanRealTC(iRoi,rwd,:) = squeeze(mean(realSubMeanTC(:,iRoi,rwd,:)));%average timecourse over subjects
        meanPermTC(iRoi,rwd,:,:) = squeeze(mean(permSubMeanTC(:,iRoi,rwd,:,:)));%average timecourse over subjects
        
        meanVar(iRoi, rwd,:) = mean(subVar(:,iRoi,rwd,:));%mean trial-by-trial variability
        stdVar(iRoi, rwd,:) = std(subVar(:,iRoi,rwd,:));%std of trial-by-trial variability
    end
    meanRealDiffStd(iRoi) = std(squeeze(meanRealTC(iRoi,1,:))) - std(squeeze(meanRealTC(iRoi,2,:)));
    for p=1:nperms
        meanPermDiffStd(iRoi,p) = std(meanPermTC(iRoi,1,p,:)) - std(meanPermTC(iRoi,2,p,:));
    end
    pVal_ffx(iRoi) = sum(meanPermDiffStd(iRoi,:) > meanRealDiffStd(iRoi))/nperms;
    
    %FFT FFX
    for rwd = 1:2
        meanRealRun(iRoi,rwd,:) = squeeze(mean(realSubMeanRun(:,iRoi,rwd,:)));%average timecourse over subjects
        meanPermRun(iRoi,rwd,:,:) = squeeze(mean(permSubMeanRun(:,iRoi,rwd,:,:)));%average timecourse over subjects
    end
    temp = abs(fft(meanRealRun(iRoi,1,:))) - abs(fft(meanRealRun(iRoi,2,:)));
    meanRealDiffFFT(iRoi) = sum(temp(ntrials+1:ntrials:(1+end/2)));
%     std(squeeze(meanRealTC(iRoi,1,:))) - std(squeeze(meanRealTC(iRoi,2,:)));
    for p=1:nperms
        temp = abs(fft(meanPermRun(iRoi,1,p,:))) - abs(fft(meanPermRun(iRoi,2,p,:)));
        meanPermDiffFFT(iRoi,p) = sum(temp(ntrials+1:ntrials:(1+end/2)));
%         meanPermDiffStd(iRoi,p) = std(meanPermTC(iRoi,1,p,:)) - std(meanPermTC(iRoi,2,p,:));
    end
    pVal_fft_ffx(iRoi) = sum(meanPermDiffFFT(iRoi,:) > meanRealDiffFFT(iRoi))/nperms;
    
    %variability
    permVar(iRoi,:,:,:) = mean(permSubVar(:,iRoi,:,:,:));%mean over subjects. permStd(iRoi,iRoi,rwd,p)
    for t=1:trialLength
       pVal_var_timecourse(iRoi,t) = sum( (permVar(iRoi,2,:,t)-permVar(iRoi,1,:,t)) > (meanVar(iRoi,2,t)-meanVar(iRoi,1,t)))/nperms;
    end
    pVal_var(iRoi) = sum( mean(permVar(iRoi,2,:,:)-permVar(iRoi,1,:,:),4) > mean(meanVar(iRoi,2,:)-meanVar(iRoi,1,:),3) )/nperms;

end
pVal_rfx
pVal_ffx
pVal_fft_rfx
pVal_fft_ffx
pVal_var_timecourse*trialLength
pVal_var
%% bin voxels within each ROI according to eccentricity, then average within bins
%bins by log eccentricity
eccMin = 0.2;
eccMax = 70;
nbins = 5;
% binBorders = (linspace(eccMin^3,eccMax^3,nbins+1))^(1/3)
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
binBorders = [0 1.3 3 7 20 100];
binBorders = [0 0.4 1.3 2 3 5 7 12 20 40 100];
nbins = length(binBorders)-1;
clear subBinTrialResponse subBinResponse binVoxels numBinVoxels meanSubVoxCorr meanVoxCorr



%%
histbins=1000;
rows=length(roiNames);
cols = length(goodSubs);
i=i+1;
figure(i)
clf
for iSub = 1:length(goodSubs)
    for iRoi= ROIs
        subplot(rows,cols,iSub + (iRoi-1)*cols)
        histogram(squeeze(permSubDiff(iSub,iRoi,:)),histbins); hold all
        vline(prctile(permSubDiff(iSub,iRoi,:),95),'k');
        vline(realSubDiff(iSub,iRoi),'r');
    end
end
        
        
%%
i=i+1;
figure(i)
clf
rows=length(roiNames);
cols = 2;
for iRoi= ROIs
    %RFX
    subplot(rows,cols,1+(iRoi-1)*cols)
        histogram(squeeze(meanPermDiff(iRoi,:)),histbins); hold all
        vline(prctile(meanPermDiff(iRoi,:),95),'k');
        vline(meanRealDiff(iRoi),'r');
        title(['random effect: ' roiNames{iRoi}]);
    %FFX
        subplot(rows,cols,2+(iRoi-1)*cols)
        histogram(squeeze(meanPermDiffStd(iRoi,:)),histbins); hold all
        vline(prctile(meanPermDiffStd(iRoi,:),95),'k');
        vline(meanRealDiffStd(iRoi),'r');
        title(['fixed effect: ' roiNames{iRoi}]);
end


%% ANOVA
clear roiSubResp
for iRoi=ROIs
    for rwd=1:2
        roiSubResp(iRoi,rwd,:) = std(realSubMeanTC(:,iRoi,rwd,:));
    end
end

varNames = {'rightBensonH', 'rightBensonL', 'rightDmnH', 'rightDmnL'};
varNames = {'Y1','Y2','Y3','Y4'};
data = squeeze([roiSubResp(2,1,:) roiSubResp(2,2,:) roiSubResp(4,1,:) roiSubResp(4,2,:)])';
tbl = array2table(data, 'VariableNames',varNames);

reward = {'H','L','H','L'};
roi = {'benson', 'benson', 'dmn','dmn'};


within = table(reward',roi', 'VariableNames',{'reward','roi'});

rm = fitrm(tbl, 'Y1-Y4~1', 'WithinDesign', within);
ranovatbl = ranova(rm, 'WithinModel', 'reward*roi')


%% Correlation between amplitude difference and FFT difference
realSubVar = mean(subVar,4);
realSubVarDiff = realSubVar(:,:,2) - realSubVar(:,:,1);
[r p] = corr(realSubDiff, realSubFftDiff)
% [r p] = corr(realSubDiff, realSubVarDiff)
% [r p] = corr(realSubFftDiff, realSubVarDiff)


toc