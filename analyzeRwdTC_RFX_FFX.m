dataFolder = '/Volumes/MH02086153MACDT-Drobo/decodingAnalysis/rwd/';
onlyCorrect=0;
nperms=10000;
onlyCorrectString = '';
if onlyCorrect
    onlyCorrectString = '_correct';
end
load([dataFolder 'rwdTC' onlyCorrectString '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', 'meanResponse', 'stdResponse',...
    'roiTC', 'allTrials', ...
    'subdirs', 'roiNames','subTrialResponse','subRunResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'rwdPupil','meanPupil','expName');

clear subMeanResponse trialStd meanTrialStd runStd meanRunStd subRwdStd trialFFTamp trialFFTphase meanTrialFFTamp meanTrialFFTphase 
clear stdTrialFFTphase runFFTamp runFFTphase meanRunFFTamp meanRunFFTphase stdRunFFTphase subFFTamp subFFTphase
clear trialMinMax meanTrialMinMax runMinMax meanRunMinMax subRwdMinMax legendCellArray
clear stdRunMax stdTrialMax
clear groupLabels pVal
clear yError yMin yMax permDiff meanPermDiff
clear permSubDiff realSubDiff permSubResp realSubResp meanPermTC permSubMeanTC realSubMeanTC meanRealTC


% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
ROIs = [1 2 7 8];
ROIs = [1 2];



%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;



%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd

% subject 4 = s0022
% subject 16 = outlier of RT and RT std
% goodSubs = [1:3 5:18];%excluding subject 22.     0.0830    0.0402
goodSubs = [1:3 5:6 8:9 11 13:18];%excluding subject 22 and all double
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




%% PERMUTATIONS
tic

for iSub = 1:length(goodSubs)%length(subdirs)
    for rwd=1:2

        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!
        if rwd==1
            firstTrial(rwd)=1;
        else %rwd==2
            firstTrial(rwd)=numTrials(iSub,1)+1;
        end
    end
    for iRoi=ROIs%length(roiNames)
        roiTrials{iRoi} = [subTrialResponse{goodSubs(iSub),iRoi,1} subTrialResponse{goodSubs(iSub),iRoi,2}];
    end
    
    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));
        for iRoi= ROIs%length(roiNames)
            for rwd=1:2
                permTrials = roiTrials{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));% 10 timepoints X #trials
                permSubMeanTC(iSub,iRoi,rwd,p,:) = mean(permTrials,2);%average timecourse over trials
                permSubResp(iSub,iRoi,rwd,p) = std(squeeze(permSubMeanTC(iSub,iRoi,rwd,p,:)));%amplitude of task response

            end
        end
    end
    
    for iRoi= ROIs%length(roiNames)
        for rwd = 1:2
            realSubMeanTC(iSub,iRoi,rwd,:) = mean(subTrialResponse{goodSubs(iSub),iRoi,rwd},2);
        end
        realSubDiff(iSub,iRoi) = std(realSubMeanTC(iSub,iRoi,1,:)) - std(realSubMeanTC(iSub,iRoi,2,:));%difference between high std and low std
        permSubDiff(iSub,iRoi,:) = squeeze(permSubResp(iSub,iRoi,1,:) - permSubResp(iSub,iRoi,2,:));%difference between high std and low std
        pVal_sub(iSub, iRoi) = sum(permSubDiff(iSub,iRoi,:) > realSubDiff(iSub,iRoi))/nperms;
        

    end
end
hypoth = pVal_sub<0.05;
% pVal_sub(:,ROIs);
for iRoi= ROIs
    %RFX
    meanPermDiff(iRoi,:) = mean(permSubDiff(:,iRoi,:));%average over subjects
    meanRealDiff(iRoi) = mean(realSubDiff(:,iRoi));%average over subjects
    pVal_rfx(iRoi) = sum(meanPermDiff(iRoi,:) > meanRealDiff(iRoi))/nperms;
    
    %FFX
    for rwd = 1:2
        meanRealTC(iRoi,rwd,:) = squeeze(mean(realSubMeanTC(:,iRoi,rwd,:)));%average timecourse over subjects
        meanPermTC(iRoi,rwd,:,:) = squeeze(mean(permSubMeanTC(:,iRoi,rwd,:,:)));%average timecourse over subjects
    end
    meanRealDiffStd(iRoi) = std(squeeze(meanRealTC(iRoi,1,:))) - std(squeeze(meanRealTC(iRoi,2,:)));
    for p=1:nperms
        meanPermDiffStd(iRoi,p) = std(meanPermTC(iRoi,1,p,:)) - std(meanPermTC(iRoi,2,p,:));
    end
    pVl_ffx(iRoi) = sum(meanPermDiffStd(iRoi,:) > meanRealDiffStd(iRoi))/nperms;
end
pVal_rfx
pVl_ffx

%%
nbins=1000;
rows=length(ROIs);
cols = length(goodSubs);
figure(1)
clf
for iSub = 1:length(goodSubs)
    for iRoi= ROIs
        subplot(rows,cols,iSub + (iRoi-1)*cols)
        hist(squeeze(permSubDiff(iSub,iRoi,:)),nbins); hold all
        vline(prctile(permSubDiff(iSub,iRoi,:),95),'k');
        vline(realSubDiff(iSub,iRoi),'r');
    end
end
        
        

figure(2)
clf
rows=length(ROIs);
cols = 2;
for iRoi= ROIs
    %RFX
    subplot(rows,cols,1+(iRoi-1)*cols)
        hist(squeeze(meanPermDiff(iRoi,:)),nbins); hold all
        vline(prctile(meanPermDiff(iRoi,:),95),'k');
        vline(meanRealDiff(iRoi),'r');
        title(['random effect: ' roiNames{iRoi}]);
    %FFX
        subplot(rows,cols,2+(iRoi-1)*cols)
        hist(squeeze(meanPermDiffStd(iRoi,:)),nbins); hold all
        vline(prctile(meanPermDiffStd(iRoi,:),95),'k');
        vline(meanRealDiffStd(iRoi),'r');
        title(['fixed effect: ' roiNames{iRoi}]);
end
%%

toc
