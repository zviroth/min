close all; clear all;
saveFolder = '/Volumes/MH02086153MACDT-Drobo/allMinBehavioral/';
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
end
tic
load([saveFolder 'behavioralData' onlyCorrectString '.mat'], 'subFolders', 'subPupil', 'runSize', 'rwdPupil',...
    'meanPupil','stdRwd','diffPupil','expName',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh','subRTvar','rwdLevel','numTrials',...
    'trialCorrectness' ,'trialResponse','trialRT','propCorrect','stairThresh');
%% mean + STD of thresholds
subThresh = mean(subMeanThresh,2);%mean over rwd, which is already averaged over staircases and runs
mean(subThresh)
std(subThresh)
keyboard
%%

plotColors = { [0 0 1],[1 0 0],[0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
goodSubs = 1:length(subFolders);
minLength=3000;
nperms=1000;
baseT = 50;
%%
clear permSubMeanTC
for iSub = 1:length(subFolders)%length(subdirs)
    for rwd=1:2
        numTrials(iSub,rwd) = size(rwdPupil{iSub,rwd},1);%the saved numTrials is for all trials, including those with no response
    end
    
    trialsPerRun(iSub) = size(subPupil{iSub,1},1);
    firstTrial(1)=1;
    firstTrial(2)=numTrials(iSub,1)+1;
    subTrials = [rwdPupil{iSub,1}(:,1:minLength); rwdPupil{iSub,2}(:,1:minLength)]; 
    subRT = [trialRT{iSub,1}; trialRT{iSub,2}]; 
    
    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));
            for rwd=1:2
                permTrials = subTrials(randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1),:);% #trials X timepoints
                permSubMeanTC(iSub,rwd,p,:) = nanmean(permTrials(:,1:minLength));
                permSubMeanStd(iSub,rwd,p) = nanstd(squeeze(permSubMeanTC(iSub,rwd,p,:)));
                permSubMeanBaseline(iSub,rwd,p) = nanmean(squeeze(permSubMeanTC(iSub,rwd,p,1:baseT)));
                %phasic and tonic for each trial and permutation
                permSubBaseline{iSub,rwd,p} = nanmean(permTrials(:,1:baseT),2);
                permSubStd{iSub,rwd,p} = nanstd(permTrials,0,2);
                permSubCorr(iSub,rwd,p) = corr(permSubStd{iSub,rwd,p},permSubBaseline{iSub,rwd,p},'rows','complete');

                permSubVarTC = nanstd(permTrials(:,1:minLength));
                permSubVar(iSub,rwd,p) = nanmean(permSubVarTC);
                
                %RT variability
                permRT = subRT(randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1),:);
                permSubRTvar(iSub,rwd,p) = std(permRT);
                
                
                keyboard
                %FFT
                f = fft(permTrials,[],2);
                ph = angle(f(:,2));
                permSubMeanPh(iSub,rwd,p) = circ_mean(ph);
                permSubPhVar(iSub,rwd,p) = circ_std(ph);
                
                f=fft(squeeze(permSubMeanTC(iSub,rwd,p,:)));
                ph = angle(f(2));
            end
    end
    for rwd=1:2
       subRwdTrials = subTrials(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1,1:minLength);
       realSubMeanTC(iSub,rwd,:) =  nanmean(subRwdTrials);
       realSubStd(iSub,rwd) = nanstd(squeeze(realSubMeanTC(iSub,rwd,:)));
       realSubBaseline(iSub,rwd) = nanmean(squeeze(realSubMeanTC(iSub,rwd,1:baseT)));
       realSubTrialStd{iSub,rwd} = nanstd(subRwdTrials,0,2);
       realSubTrialBaseline{iSub,rwd} = nanmean(subRwdTrials(:,1:baseT),2);
%        realSubRunStd{iSub,rwd} = mean(reshape(realSubTrialStd{iSub,rwd},trialsPerRun(iSub),[]),2);%mean across runs, phasic timeline
%        realSubRunBaseline{iSub,rwd} = mean(reshape(realSubTrialBaseline{iSub,rwd},trialsPerRun(iSub),[]),2);%mean across runs, tonic timeline
       realSubCorr(iSub,rwd) = corr(squeeze(realSubTrialStd{iSub,rwd}),squeeze(realSubTrialBaseline{iSub,rwd}),'rows','complete');
       realSubMeanStd(iSub,rwd) = mean(realSubTrialStd{iSub,rwd});
       realSubMeanBaseline(iSub,rwd) = mean(realSubTrialBaseline{iSub,rwd});
       
       realSubVarTC = nanstd(subRwdTrials);
       realSubVar(iSub,rwd) = nanmean(realSubVarTC);
       
       realSubRTvar(iSub,rwd) = std(trialRT{iSub,rwd});
       
       pval_corr_sub(iSub,rwd) = sum(permSubCorr(iSub,rwd,:) > realSubCorr(iSub,rwd))/nperms;
    end
    realSubStdDiff(iSub) = realSubStd(iSub,1) - realSubStd(iSub,2);
    realSubBaselineDiff(iSub) = realSubBaseline(iSub,1) - realSubBaseline(iSub,2);
    permSubStdDiff(iSub,:) = permSubMeanStd(iSub,1,:) - permSubMeanStd(iSub,2,:);
    permSubBaselineDiff(iSub,:) = permSubMeanBaseline(iSub,1,:) - permSubMeanBaseline(iSub,2,:);
    
    %     pVal_sub(iSub, iRoi) = sum(permSubDiff(iSub,iRoi,:) > realSubDiff(iSub,iRoi))/nperms;

end
%%
%average tonic and phasic timelines across subjects with same # of trials
%42 trials: 1:10, 16 trials: 11, 25 trials: 12:13

% runSubList = 1:10;
% clear subRunPhasic subRunTonic
% for iSub=runSubList
%    for rwd=1:2
%       subRunPhasic(iSub,rwd,:) = realSubRunStd{iSub,rwd};
%       subRunTonic(iSub,rwd,:) = realSubRunBaseline{iSub,rwd};
%    end
% end
% runPhasic = squeeze(nanmean(subRunPhasic));
% runTonic = squeeze(nanmean(subRunTonic));
% 
% figure(1)
% clf
% subplot(1,2,1)
% plot(zscore(runPhasic,0,2)')
% hold on
% subplot(1,2,2)
% plot(zscore(runTonic,0,2)')


%%
meanPermStdDiff = nanmean(permSubStdDiff);%average over subjects
meanRealStdDiff = nanmean(realSubStdDiff);
meanPermBaselineDiff =nanmean(permSubBaselineDiff);%average over subjects
meanRealBaselineDiff = nanmean(realSubBaselineDiff);

meanRealVar = nanmean(realSubVar);
meanRealVarDiff = meanRealVar(2)-meanRealVar(1);
meanPermVar = squeeze(nanmean(permSubVar));
meanPermVarDiff = meanPermVar(2,:) - meanPermVar(1,:);
pVal_var = sum(meanPermVarDiff > meanRealVarDiff)/nperms;

%RFX
pVal_std_rfx = sum(meanPermStdDiff >= meanRealStdDiff)/nperms;
pVal_baseline_rfx = sum(meanPermBaselineDiff >= meanRealBaselineDiff)/nperms;
%FFX
for rwd = 1:2
    meanRealTC(rwd,:) = squeeze(nanmean(realSubMeanTC(:,rwd,:)));%average timecourse over subjects
    meanPermTC(rwd,:,:) = squeeze(nanmean(permSubMeanTC(:,rwd,:,:)));%average timecourse over subjects
end
realStdDiff = nanstd(squeeze(meanRealTC(1,:))) - nanstd(squeeze(meanRealTC(2,:)));
realBaselineDiff = nanmean(squeeze(meanRealTC(1,1:baseT))) - nanmean(squeeze(meanRealTC(2,1:baseT)));

for p=1:nperms
    permStdDiff(p) = nanstd(meanPermTC(1,p,:)) - nanstd(meanPermTC(2,p,:));
    permBaselineDiff(p) = nanmean(meanPermTC(1,p,1:baseT)) - nanmean(meanPermTC(2,p,1:baseT));
end
pVal_std_ffx = sum(permStdDiff >= realStdDiff)/nperms;
pVal_baseline_ffx = sum(permBaselineDiff >= realBaselineDiff)/nperms;
    
    

pVal_std_rfx
pVal_baseline_rfx
pVal_std_ffx
pVal_baseline_ffx

pVal_var



permRTvar = squeeze(mean(permSubRTvar));
permRTvarDiff = permRTvar(2,:)-permRTvar(1,:);
realRTvar = mean(realSubRTvar);
realRTvarDiff = realRTvar(2) - realRTvar(1);
pVal_RTvar = sum(permRTvarDiff >= realRTvarDiff)/nperms;

pVal_RTvar
%%

% goodSubs = [1:4 6:13];
goodSubs = 1:13;
for rwd=1:2
   [realPhasicTonicCorr(rwd) corrPval(rwd)]= corr( realSubBaseline(goodSubs,rwd),realSubStd(goodSubs,rwd));
end
% for p=1:nperms
%     for rwd=1:2
%         permPhasicTonicCorr(rwd,p)= corr(permSubMeanBaseline(goodSubs,rwd,p),permSubMeanStd(goodSubs,rwd,p)); 
%     end
% end
% for rwd=1:2
%     pVal_corr(rwd) = sum(permPhasicTonicCorr(rwd,:) > realPhasicTonicCorr(rwd))/nperms;
%     pVal_corrInverse(rwd) = sum(permPhasicTonicCorr(rwd,:) < realPhasicTonicCorr(rwd))/nperms;
% end
% 
% pVal_corr
% pVal_corrInverse
toc
%%
save([saveFolder 'RTvariability' onlyCorrectString '.mat'], 'realRTvarDiff','permRTvarDiff','pVal_RTvar',...
    'realRTvar','realSubRTvar','trialRT','permRTvar','permSubRTvar','subFolders','nperms',...
    'pVal_std_rfx','pVal_baseline_rfx','pVal_var','meanRealStdDiff','realSubStdDiff','meanPermVarDiff',...
    'permBaselineDiff','permSubVar','permSubStdDiff','realBaselineDiff','realStdDiff');


%%

%%
% figure(1)
% clf
% rows=2;
% cols=2;
% for rwd=1:2
%     subplot(rows,cols,rwd)
%     plot(realSubBaseline(goodSubs,rwd),realSubStd(goodSubs,rwd),'.','markersize',20);
%     
%     subplot(rows,cols,rwd+cols)
%     histogram(permPhasicTonicCorr(rwd,:)); hold on;
%     vline(prctile(permPhasicTonicCorr(rwd,:),95),'k');
%     vline(prctile(permPhasicTonicCorr(rwd,:),5),'k');
%     vline(realPhasicTonicCorr(rwd),'r');
% end

%% ANOVA
% varNames = {'baseline', 'std'};
varNames = {'baselineH', 'baselineL','stdH','stdL'};
varNames = {'Y1','Y2','Y3','Y4'};
% varNames = {'Y1','Y2'};
% tbl = table(realSubBaseline(:,1),realSubBaseline(:,2),realSubStd(:,1),realSubStd(:,2),'VariableNames',varNames);
tbl = array2table([realSubBaseline realSubStd], 'VariableNames',varNames);
reward = {'H','L','H','L'};
pupilMeasure = {'tonic','tonic','phasic','phasic'};

within = table(reward',pupilMeasure', 'VariableNames',{'reward','pupilMeasure'});

rm = fitrm(tbl, 'Y1-Y4~1', 'WithinDesign', within);
ranovatbl = ranova(rm, 'WithinModel', 'reward*pupilMeasure')
