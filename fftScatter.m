close all
clear all
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
regressGlobalMean = 0;
ConcatProj = 0;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
TR=1.5;
nperms=100;
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
elseif onlyCorrect==0
    onlyCorrectString = '_validresponse';
end
fmriUnits = '% change image intensity';
zScoreString = '';
if toZscore
    zScoreString = '_zscored';
    fmriUnits = 'std image intensity';
end

globalMeanString = '';
if regressGlobalMean==1
    globalMeanString = '_globalRegressed';
elseif regressGlobalMean==2
    globalMeanString = '_bensonRegressed';
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
ROIs = [1 2 9 10];
ROIs = [1:6];
dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;


ntrials=15;

clear allTrials meanRwdStd binCenters

goodSubs = 1:length(subFolders);
% goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% subjects. 0.0457  0.0140,     0.0553    0.0214

for iSub = goodSubs%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(ROIs)
        for rwd=1:2
            reshapedTrials{rwd} = reshape(subTrialResponse{iSub,iRoi,rwd},10,[]);
        end
        subTrials = [reshapedTrials{1} reshapedTrials{2}];
        subMeanResponse(iSub,iRoi,:) = mean(subTrials,2);
     end
end



%%
iRoi=2;
for iSub=1:length(subFolders)
    
    f(iSub,:)=fft(subMeanResponse(iSub,iRoi,:));
    for rwd=1:2
        fRwd(iSub,rwd,:)=fft(subResponse(iSub,iRoi,rwd,:));
    end
end
subReal(:,1) = real(f(:,2));
subImag(:,1) = imag(f(:,2));

for rwd=1:2
    subReal(:,rwd+1) = real(fRwd(:,rwd,2));
    subImag(:,rwd+1) = imag(fRwd(:,rwd,2));
end

%%

i=1; figure(i); clf
titleStr = {'mean','H','L'};
% plot(subReal, subImag, 'k.');
rows=1;cols=3;
nSTD = tinv(0.99,length(subFolders)-1)/sqrt(length(subFolders));
nSTD = 2.5;
% nSTD = tinv(0.98,length(subFolders)-1);
% nSTD = 2.576/sqrt(length(subFolders));
for j=1:3
    subplot(rows,cols,j)
    scatter(subReal(:,j), subImag(:,j), 30, 'black','filled');
    hold on
    drawConfidenceEllipse([mean(subReal(:,j)) mean(subImag(:,j))],cov([subReal(:,j) subImag(:,j)]),'alpha=0.3','nSTD',nSTD);
    hline(0)
    vline(0)
    xlabel([ 'real Fourier component (' fmriUnits ')']);
    if j==1; ylabel([ 'imaginary Fourier component (' fmriUnits ')']); end
    title(titleStr{j});
    drawPublishAxis('xLabelOffset', -7/64,'yLabelOffset', -10/64, 'xAxisMargin', 2/64, 'yAxisMargin', 2/64,'xAxisMinMaxSetByTicks',0,...
        'labelFontSize',7);
    axis square
end
set(gcf,'position',[10 10 25 6]);

print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig4_fourierScatterPlot_' ConcatProjStr '.pdf']);

%%

i=i+1; figure(i); clf
titleStr = {'mean','H','L'};
% plot(subReal, subImag, 'k.');
rows=1;cols=3;
nSTD = tinv(0.99,length(subFolders)-1)/sqrt(length(subFolders));
nSTD = 2.5;
% nSTD = tinv(0.98,length(subFolders)-1);
% nSTD = 2.576/sqrt(length(subFolders));
for j=2:3
%     subplot(rows,cols,j)
    scatter(subReal(:,j), subImag(:,j), 30, plotColors{j-1},'filled');
    hold on

%     drawConfidenceEllipse([mean(subReal(:,j)) mean(subImag(:,j))],cov([subReal(:,j) subImag(:,j)]),'alpha=0.3','nSTD',nSTD);
    hline(0)
    vline(0)
    xlabel([ 'real Fourier component (' fmriUnits ')']);
    ylabel([ 'imaginary Fourier component (' fmriUnits ')']);
%     title(titleStr{j});
    drawPublishAxis('xLabelOffset', -7/64,'yLabelOffset', -10/64, 'xAxisMargin', 2/64, 'yAxisMargin', 2/64,'xAxisMinMaxSetByTicks',0,...
        'labelFontSize',7);
    axis square
end
for iSub=1:length(subFolders)
    plot( [subReal(iSub,2) subReal(iSub,3)], [subImag(iSub,2) subImag(iSub,3)],'-','color',dsSurfaceContrast*[1 1 1]);
end
    
set(gcf,'position',[10 10 25 6]);

print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig4_fourierScatterPlot2_' ConcatProjStr '.pdf']);

%%
i=i+1; figure(i); clf
chosenROIs = [2 4 6];
rows=1; cols=length(chosenROIs);
    j=1;%only mean across H and L
for r = 1:length(chosenROIs)
    iRoi = chosenROIs(r);

    for iSub=1:length(subFolders)
        f(iSub,:)=fft(subMeanResponse(iSub,iRoi,:));
        for rwd=1:2
            fRwd(iSub,rwd,:)=fft(subResponse(iSub,iRoi,rwd,:));
        end
    end
    subReal(:,1) = real(f(:,2));
    subImag(:,1) = imag(f(:,2));
    for rwd=1:2
        subReal(:,rwd+1) = real(fRwd(:,rwd,2));
        subImag(:,rwd+1) = imag(fRwd(:,rwd,2));
    end
    
    subplot(rows,cols,r)
    
    scatter(subReal(:,j), subImag(:,j), 30, 'black','filled');
    hold on
    drawConfidenceEllipse([mean(subReal(:,j)) mean(subImag(:,j))],cov([subReal(:,j) subImag(:,j)]),'alpha=0.3','nSTD',nSTD);
    hline(0)
    vline(0)
    xlabel({'real Fourier component';['(' fmriUnits ')']});
%     if r==1
        ylabel({ 'imaginary Fourier component';['(' fmriUnits ')']}); 
%     end
    title(roiNames{iRoi});
    drawPublishAxis('xLabelOffset', -7/64,'yLabelOffset', -10/64, 'xAxisMargin', 2/64, 'yAxisMargin', 2/64,'xAxisMinMaxSetByTicks',0,...
        'labelFontSize',7);
    axis square
end
set(gcf,'position',[10 10 25 6]);

print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig4_fourierScatterPlot_bothROIs_' ConcatProjStr '.pdf']);
