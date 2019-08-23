onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
regressGlobalMean = 0;
ConcatProj = 1;
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
plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotColors = { [1 0 0],[0 0 1],[0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
ROIs = [1 2 9 10];
ROIs = [1 2];



%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;



%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd binCenters

% subject 4 = s0022
% subject 16 = outlier of RT and RT std
% goodSubs = [1:3 5:18];%excluding subject 22.     0.0830    0.0402
% goodSubs = [1:3 5:6 8:9 11 13:18];%excluding subject 22 and all double
goodSubs = 1:length(subFolders);
% goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% subjects. 0.0457  0.0140,     0.0553    0.0214

for iSub = goodSubs%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(ROIs)
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
for rwd=1:2
    for iRoi = 1:length(ROIs)
        meanResponse(rwd,iRoi,:) = mean(subResponse(goodSubs,iRoi,rwd,:));
        stdResponse(rwd,iRoi,:) = std(subResponse(goodSubs,iRoi,rwd,:));
    end
end



%% PERMUTATIONS
tic


%% bin voxels within each ROI according to eccentricity, then average within bins
%bins by log eccentricity
eccMin = 0.2;
eccMax = 40;
nbins = 10;
% binBorders = (linspace(eccMin^3,eccMax^3,nbins+1))^(1/3)
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);

% binBorders = 1:12;
% binBorders = [0 1.3 3 7 20 100];
% binBorders = [0 0.4 1.3 2 3 5 7 12 20 40 100];
% binSpacing = 0.2:0.4:7;
% binBorders(1)=0;
% for i=1:length(binSpacing)
%     binBorders(i+1) = binBorders(i)+binSpacing(i);
% end
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
% binBorders=exp(linspace(log(eccMin),log(eccMax),nbins+1));
% binBorders = (linspace(eccMin,eccMax,nbins+1));
% minBinVoxels=5;%minimum voxels per bin, otherwise subject is excluded from average
clear subBinTrialResponse subBinResponse binVoxels numBinVoxels meanSubVoxCorr meanVoxCorr
for iSub = 1:length(goodSubs)%length(subdirs)
    for iRoi= ROIs
        for rwd=1:2
            temp = trialCorrectness{iSub,rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);
            temp = trialResponse{iSub,rwd}(:,2:end-1);
            trialResponseVec = temp(:);
            if onlyCorrect ==1 %ONLY CORRECT
                    goodTrials = trialCorrectnessVec==1;
                    %                 numTrials = sum(trialCorrectnessVec);
                elseif onlyCorrect ==2 % ONLY INCORRECT!!!
                    goodTrials = trialCorrectnessVec==0 & trialResponseVec>0;
                else % including all trials with a response
                    goodTrials = trialResponseVec>0;
            end
                
            
        end
    end
end

%%
dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;
i=1;
figure(i)
clf
rows=2;
cols=length(ROIs);
iSub=2;
for iRoi = ROIs
    
    subplot(rows,cols,iRoi);
    
    for rwd=1:2
        dsErrorsurface(TR*(0:trialLength-1), squeeze(subResponse(iSub,iRoi,rwd,:)), squeeze(subStd(iSub,iRoi,rwd,:))./sqrt(size(subTrialResponse{iSub,iRoi,rwd},2)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(TR*(0:trialLength-1),squeeze(subResponse(iSub,iRoi,rwd,:)), 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',20);
    end
%     xlim(TR*[0 1+trialLength]);
    title([subFolders{iSub}(1:5) ' ' roiNames{iRoi}]);
    xlabel('time (sec)');
    ylabel(['BOLD signal (' fmriUnits ')']);
    axis square
    box on
    
end
%group mean
for iRoi = ROIs
    subplot(rows,cols,iRoi+cols);
    
    for rwd=1:2
        
        dsErrorsurface(TR*(0:trialLength-1), squeeze(meanResponse(rwd,iRoi,:)), squeeze(stdResponse(rwd,iRoi,:))./sqrt(length(goodSubs)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(TR*(0:trialLength-1),squeeze(meanResponse(rwd,iRoi,:)), 'Color', plotColors{rwd}, 'linewidth', linewidth);
        
    end
%     xlim(TR*[0 1+trialLength]);
    title(roiNames{iRoi})
    xlabel('time (sec)');
    ylabel(['BOLD signal (' fmriUnits ')']);
    axis square
    box on
end

for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    drawPublishAxis('xLabelOffset', -9/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 4/64);
end
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig3_bold' ConcatProjStr '.pdf']);
%% mean response for all subjects
i=i+1;
figure(i)
clf
rows=4;
cols=ceil(length(subFolders)/rows);
iRoi=2;
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
   xlabel('time (sec)');
   ylabel(['BOLD signal (' fmriUnits ')']);
   axis square
   box on
%    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 2/64,'yAxisMinMaxSetByTicks',0,'xTick',[0 5 10 15]);
   drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -14/64, 'xAxisMargin', 4/64, 'yAxisMargin', 2/64,'xTick',[0 5 10 15]);
end
set(gcf,'position',[10 10 30 20]);

%% single subject
i=i+1;
figure(i)
clf
cols=length(subFolders);
rows=length(ROIs);
for iRoi = ROIs
    for iSub = 1:length(subFolders)
        subplot(rows,cols,(iRoi-1)*cols+iSub);
        for rwd=1:2
            plot(squeeze(subResponse(iSub,iRoi,rwd,:)), '.-', 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',20);
            hold on
        end
    end
end

for isubplot=1:rows*cols
   subplot(rows,cols,isubplot)
   xlabel('time (sec)');
   ylabel(['BOLD signal (' fmriUnits ')']);
   axis square
   box on
   drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -8/64, 'xAxisMargin', 4/64, 'yAxisMargin', 4/64,'xAxisMinMaxSetByTicks',0);
end


%% FFT of average run, for group average and single subject
% i=i+1;
% figure(i)
% clf
% cols=length(ROIs);
% rows=length(subFolders);
rows=length(ROIs);
cols=length(subFolders);

for iSub=1:length(subFolders)
    for iRoi = ROIs
%         subplot(rows,cols,iSub + (iRoi-1)*cols)
        for rwd=1:2
            subMeanRun(iRoi,rwd,:) = squeeze(subMeanRunTC(iSub,iRoi,rwd,:));
            subMeanRun(iRoi,rwd,:) = subMeanRun(iRoi,rwd,:) - mean(subMeanRun(iRoi,rwd,:));
            %         subMeanRun(iRoi,rwd,:) = zscore(subMeanRun(iRoi,rwd,:));
            runFFT = squeeze(abs(fft(subMeanRun(iRoi,rwd,:))));
            runFFT = runFFT(1:ceil((length(runFFT)+1)/2));
%             plot(0:length(runFFT)-1, runFFT,'.-','Color',plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);
%             hold on
        end
%         title([subFolders{iSub}(1:4) ' ' roiNames{iRoi}]);
%         axis square
%         box on
    end
end

%% FFT of average run, for group average
i=i+1;
figure(i)
clf
rows=2;
cols=length(ROIs);
meanRunTC = squeeze(mean(subMeanRunTC));
ymax=0;
% FFT of average run for s0008

iSub=2;

for iRoi = ROIs
    subplot(rows,cols,iRoi)
    for rwd=1:2
        subMeanRun(iRoi,rwd,:) = squeeze(subMeanRunTC(iSub,iRoi,rwd,:));
        subMeanRun(iRoi,rwd,:) = subMeanRun(iRoi,rwd,:) - mean(subMeanRun(iRoi,rwd,:));
        %         subMeanRun(iRoi,rwd,:) = zscore(subMeanRun(iRoi,rwd,:));
        runFFT = squeeze(abs(fft(subMeanRun(iRoi,rwd,:))));
        % Compute mean fourier amplitude

        
        nframes = length(runFFT);
        %from plotMeanFourierAmp.m
        runFFT = runFFT / (nframes/2);
        frequencies = [0:nframes-1]/(nframes*TR);
        runFFT = runFFT(1:floor(nframes/2));
        frequencies = frequencies(1:floor(nframes/2));
        


%         runFFT = runFFT./TR;%is this correct?
%         runFFT = runFFT(1:ceil((length(runFFT)+1)/2));
%         plot(0:length(runFFT)-1, runFFT,'.-','Color',plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);
        plot(frequencies, runFFT,'.-','Color',plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);

        hold on
            if max(runFFT(:))>ymax
                ymax= max(runFFT(:));
            end
        %         plot(1:size(subMeanRun,3), squeeze(subMeanRun(iRoi,rwd,:)),  '.-', 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);
    end
    title([subFolders{iSub}(1:4) ' ' roiNames{iRoi}]);
    axis tight
    axis square
    box on
end
%group mean
for iRoi = ROIs
        subplot(rows,cols,iRoi+cols)
        for rwd=1:2
            meanRunTC(iRoi,rwd,:) = meanRunTC(iRoi,rwd,:) - mean(meanRunTC(iRoi,rwd,:));
            runFFT = squeeze(abs(fft(meanRunTC(iRoi,rwd,:))));
            
            nframes = length(runFFT);
            %from plotMeanFourierAmp.m
            runFFT = runFFT / (nframes/2);
            frequencies = [0:nframes-1]/(nframes*TR);
            runFFT = runFFT(1:floor(nframes/2));
            frequencies = frequencies(1:floor(nframes/2));
        
            
%         plot(0:length(runFFT)-1, runFFT,'.-','Color',plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);
            plot(frequencies, runFFT,'.-','Color',plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);
            hold on
            if max(runFFT(:))>ymax
                ymax= max(runFFT(:));
            end
        end
        title(roiNames{iRoi});
        axis square
        box on
end

%%
for isubplot=1:rows*cols
   subplot(rows,cols,isubplot)
   ylim([0 ymax]);
   xlabel('frequency (cycles/run)');
   ylabel(['Fourier amplitude  (' fmriUnits ')']);
   drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 2/64, 'yAxisMargin', 2/64,'xAxisMinMaxSetByTicks',0);
end
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig3_fft' ConcatProjStr '.pdf']);
%%
toc



