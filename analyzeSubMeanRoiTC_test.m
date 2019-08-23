dataFolder = '/Volumes/MH02086153MACDT-Drobo/decodingAnalysis/rwd/';
onlyCorrect=1;
onlyCorrectString = '';
if onlyCorrect
    onlyCorrectString = '_correct';
end
load([dataFolder 'subMeanRoiTC' onlyCorrectString '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', 'meanResponse', 'stdResponse',...
    'roiTC', 'allTrials', ...
    'subdirs', 'roiNames','subTrialResponse','subRunResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'rwdPupil','meanPupil','expName');

clear subMeanResponse trialStd meanTrialStd runStd meanRunStd subRwdStd trialFFTamp trialFFTphase meanTrialFFTamp meanTrialFFTphase 
clear stdTrialFFTphase runFFTamp runFFTphase meanRunFFTamp meanRunFFTphase stdRunFFTphase subFFTamp subFFTphase
clear trialMinMax meanTrialMinMax runMinMax meanRunMinMax subRwdMinMax legendCellArray
clear stdRunMax stdTrialMax
clear groupLabels pVal
clear yError yMin yMax


% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
ROIs = [1 2 7 8];
ROIs = [1 2];

% ROIs = 1:4;

for i=1:length(ROIs)
    groupLabels{i} = roiNames{ROIs(i)};
end
legendHighLow = [];
for iRoi= 1:length(groupLabels)
    legendHighLow{(iRoi-1)*2+1} = [groupLabels{iRoi} ': H'];
    legendHighLow{iRoi*2} = [groupLabels{iRoi} ': L'];
end
% ROIs = [1:2];

% roiNames(5:end) = [];
% roiNames([1:2 4:end]) = [];
% roiNames(end) = [];
%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;
plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
figure(1)
clf
cols=ceil(length(subdirs)/2);
rows = ceil(length(subdirs)/cols);

for iSub = 1:length(subdirs)
    subplot(rows,cols,iSub)
    for iRoi= ROIs%1:length(roiNames)
        
        for rwd=1:2
            
%             STD
            subRwdStd(iSub,iRoi,rwd) = std(squeeze(subResponse(iSub,iRoi,rwd,:)));

            % plot mean timecourse
            plot(squeeze(subResponse(iSub,iRoi,rwd,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
            title([subdirs(iSub).name(1:5) newline expName{iSub,rwd}]);
            hold on
        end
        
    end

end
legend(legendHighLow);
set(gcf,'position',[250 500 1200 600]);


%% average timecourse for all ROIs
figure(5)
clf
cols=ceil(length(subdirs)/2);
rows = ceil(length(subdirs)/cols);
for iRoi= ROIs%1:length(roiNames)
    for rwd=1:2
        plot(squeeze(meanResponse(rwd,iRoi,:)),plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
        hold on
    end
end
legend(legendHighLow);



%% PERMUTATIONS
% first combine all trials across all subjects
clear allTrials
goodSubs = [1:3 5:18];
goodSubs = 1:18;
for iSub = goodSubs%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(roiNames)
        for rwd=1:2
            reshapedTrials = reshape(subTrialResponse{iSub,iRoi,rwd},10,[]);
            if iSub==1
                allTrials{iRoi,rwd} = reshapedTrials;
            else
                allTrials{iRoi,rwd} = [allTrials{iRoi,rwd} reshapedTrials];
            end
        end
     end
end

nperms=10000;
for rwd=1:2
    numTrials(rwd) = size(allTrials{1,rwd},2);
    if rwd==1
        firstTrial(rwd)=1;
    else %rwd==2
        firstTrial(rwd)=numTrials(1)+1;
    end
end
for iRoi=ROIs%length(roiNames)
   roiTrials{iRoi} = [allTrials{iRoi,1} allTrials{iRoi,2}]; 
end

for p=1:nperms
    randOrder = randperm(numTrials(1)+numTrials(2));
    for iRoi= ROIs%length(roiNames)
        for rwd=1:2
            permTrials = roiTrials{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(rwd)-1));
            permStd(iRoi,rwd,p) = std(mean(permTrials,2));
        end
    end
end
for iRoi= ROIs%length(roiNames)
    realDiff(iRoi) = std(mean(allTrials{iRoi,rwd},1)) - std(mean(allTrials{iRoi,rwd},2)); 
%     realDiff(iRoi) = meanRwdStd(iRoi,1) - meanRwdStd(iRoi,2);
%     realDiff = stdAmplitude(iRoi,1) - stdAmplitude(iRoi,2);
    permDiff = squeeze(permStd(iRoi,1,:) - permStd(iRoi,2,:));
    pVal(iRoi) = sum(permDiff>realDiff(iRoi))/nperms;
end
    
pVal(ROIs)

for iRoi = ROIs
   [h pVal_ttest(iRoi)] = ttest(subRwdStd(:,iRoi,1),subRwdStd(:,iRoi,2));
end
%%
for iSub = 1:length(subdirs)
   for rwd=1:2
       meanPropCorrect(iSub,rwd) = mean(propCorrect{iSub,rwd});
   end
end

