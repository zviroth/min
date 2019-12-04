close all
clear all
tic
onlyCorrect=0;%1=correct,2=incorrect,0=all
toZscore=1;%0 or 1
regressGlobalMean = 0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';

%%
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
%%
load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'globalMean','regressBetasGlobal','runRwd',...
    'subRoiRuns','runMeanFFT',...
    'allVoxTrialResponse','allVoxTaskPhase','allVoxTaskAmp','allVoxTaskCo',...
    'voxTrials','voxGoodTrials','meanVoxTrial');

numSubs=length(subFolders);

for iSub=1:numSubs
    for iRoi=1:length(roiNames)
        for rwd=1:2

            meanVoxTrial{iSub,iRoi,rwd};%vox,t
            stdAmp{iSub,iRoi}(rwd,:) = squeeze(std(meanVoxTrial{iSub,iRoi,rwd},0,2));
            temp = fft(meanVoxTrial{iSub,iRoi,rwd},[],2);
            ph{iSub,iRoi}(rwd,:) = angle(temp(:,2));
            amp{iSub,iRoi}(rwd,:) = abs(temp(:,2));
            
            voxGoodTrials{iSub,iRoi,rwd}; %vox,t,trial
            stdAmpTrial{iSub,iRoi,rwd} = squeeze(std(voxGoodTrials{iSub,iRoi,rwd},0,2));%vox,trial
            temp = fft(voxGoodTrials{iSub,iRoi,rwd},[],2);
            phTrial{iSub,iRoi,rwd} = squeeze(angle(temp(:,2,:)));
            ampTrial{iSub,iRoi,rwd} = squeeze(abs(temp(:,2,:)));
            
            phVar{iSub,iRoi}(rwd,:) = squeeze(circ_std(phTrial{iSub,iRoi,rwd},[],[],2));
            stdAmpVar{iSub,iRoi}(rwd,:) = squeeze(nanstd(stdAmpTrial{iSub,iRoi,rwd},0,2));
            ampVar{iSub,iRoi}(rwd,:) = squeeze(nanstd(ampTrial{iSub,iRoi,rwd},0,2));
        end
                    
        stdAmpDiff{iSub,iRoi} = stdAmp{iSub,iRoi}(1,:) - stdAmp{iSub,iRoi}(2,:);
        phDiff{iSub,iRoi} = circ_dist(ph{iSub,iRoi}(1,:),ph{iSub,iRoi}(2,:));
        ampDiff{iSub,iRoi} = amp{iSub,iRoi}(1,:) - amp{iSub,iRoi}(2,:);
        
        stdAmpVarDiff{iSub,iRoi} = stdAmpVar{iSub,iRoi}(1,:) - stdAmpVar{iSub,iRoi}(2,:);
        phVarDiff{iSub,iRoi} = phVar{iSub,iRoi}(1,:) - phVar{iSub,iRoi}(2,:);        
        ampVarDiff{iSub,iRoi} = ampVar{iSub,iRoi}(1,:) - ampVar{iSub,iRoi}(2,:);      
       
        concatForCorr = [stdAmpDiff{iSub,iRoi};  phDiff{iSub,iRoi}; ampDiff{iSub,iRoi}; ...
            stdAmpVarDiff{iSub,iRoi}; phVarDiff{iSub,iRoi}; ampVarDiff{iSub,iRoi}];
        concatForCorr = concatForCorr(:,~isnan(phDiff{iSub,iRoi}));
        
        [pearsonsCorr(iSub,iRoi,:,:), pearsonsPval(iSub,iRoi,:,:)] = corr(concatForCorr');

    end    
end
%%
labels = {'\Deltastd','\Deltaph','\Deltaamp','\Deltavar(std)','\Deltavar(ph)','\Deltavar(amp)'};
figure(1); clf
rows=2;
cols=2;
meanPearsonsCorr= squeeze(nanmean(pearsonsCorr));%mean over subjects
for iRoi=1:length(roiNames)
   subplot(rows,cols,iRoi)
   imagesc(squeeze(meanPearsonsCorr(iRoi,:,:)));
   xticklabels(labels);
   yticklabels(labels);
end

figure(2); clf
meanPearsonsPval= squeeze(nanmean(pearsonsPval));%mean over subjects
for iRoi=1:length(roiNames)
   subplot(rows,cols,iRoi)
   imagesc(squeeze(meanPearsonsPval(iRoi,:,:)));
   xticklabels(labels);
   yticklabels(labels);
end

%% plot polar histogram of task related phase
figure(2); clf
rows=2;
cols=numSubs;
nbins=20;
nVox=1000;
thresh=0.015;
for iSub=1:numSubs
%     bestVox = sort(allVoxTaskCo{iSub,1}+allVoxTaskCo{iSub,2},'descend');
    for rwd=1:2
        subplot(rows,cols,iSub+(rwd-1)*cols)
%         polarhistogram(allVoxTaskPhase{iSub,rwd}(bestVox(1:nVox)),nbins);
polarhistogram(allVoxTaskPhase{iSub,rwd}(allVoxTaskCo{iSub,rwd}>thresh),nbins);
        set(gca,'RTickLabel',[]);
        set(gca,'ThetaTickLabel',[]);
    end
end
