onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
regressGlobalMean = 0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
TR=1.5;
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
plotColors = { [1 0 0],[0 0 1],[0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
ROIs = [1 2 9 10];
ROIs = [1:4];



%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;



%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd binCenters

goodSubs = 1:length(subFolders);
% goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% subjects. 0.0457  0.0140,     0.0553    0.0214


for rwd=1:2
    for iRoi = 1:length(ROIs)
        meanResponse(rwd,iRoi,:) = mean(subResponse(goodSubs,iRoi,rwd,:));
        stdResponse(rwd,iRoi,:) = std(subResponse(goodSubs,iRoi,rwd,:));
    end
end

% minBinVoxels=5;%minimum voxels per bin, otherwise subject is excluded from average
clear subBinTrialResponse subBinResponse binVoxels numBinVoxels meanSubVoxCorr meanVoxCorr

%%
dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;


%% mean response for all subjects
i=1;
figure(i)
clf
rows=3;
cols=ceil(length(subFolders)/rows);
iRoi=4;
for iSub=1:length(subFolders)
    
        
        subplot(rows,cols,iSub);
        
        for rwd=1:2
            dsErrorsurface(TR*(0:trialLength-1), squeeze(subResponse(iSub,iRoi,rwd,:)), squeeze(subStd(iSub,iRoi,rwd,:))./sqrt(size(subTrialResponse{iSub,iRoi,rwd},2)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
            hold on
            plot(TR*(0:trialLength-1),squeeze(subResponse(iSub,iRoi,rwd,:)), 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',20);
        end
%         xlim(TR*[0 1+trialLength]);
        title([subFolders{iSub}(1:5) ' ' roiNames{iRoi}]);
end
for isubplot=1:rows*cols
   subplot(rows,cols,isubplot)
   title(['S' num2str(isubplot)]);
   xlabel('time (sec)');
   if mod(isubplot,cols)==1
       t = {'BOLD signal'; ['(' fmriUnits ')'] };
       ylabel(t);
   end
   axis square

%    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 2/64,'yAxisMinMaxSetByTicks',0,'xTick',[0 5 10 15]);
        drawPublishAxis('xLabelOffset', -9/64,'yLabelOffset', -18/64, 'xAxisMargin', 2/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',0,...
       'labelFontSize',7);
end
set(gcf,'position',[10 10 21 15]);


print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig4_' roiNames{iRoi} '_' ConcatProjStr '.pdf']);
% print('-painters','-dpdf','~/Documents/MATLAB/min/figures/fig4_allSubjects_rightDMN.pdf');

