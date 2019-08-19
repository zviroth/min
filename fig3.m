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
ROIs = 1:4;
ROIs = [2];

cols=9;
rows=length(ROIs);
subplots = {1:3, 4:6, 7 , 9};
ntrials=15;


clear allTrials meanRwdStd binCenters

% goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% subjects. 0.0457  0.0140,     0.0553    0.0214
figure(1)
clf
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

%% SINGLE SUBJECT AVERAGE TRIAL
dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;
i=1;
figure(i)
clf

iSub=2;
for r = 1:length(ROIs)
    iROI = ROIs(r);
%     subplot(rows,cols,(r-1)*cols+1)
    subplot(rows,cols,subplots{1})
    for rwd=1:2
        dsErrorsurface(TR*(0:trialLength-1), squeeze(subResponse(iSub,iRoi,rwd,:)), squeeze(subStd(iSub,iRoi,rwd,:))./sqrt(size(subTrialResponse{iSub,iRoi,rwd},2)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(TR*(0:trialLength-1),squeeze(subResponse(iSub,iRoi,rwd,:)), 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',20);
    end
%     xlim(TR*[0 1+trialLength]);
%     title([subFolders{iSub}(1:5) ' ' roiNames{iRoi}]);
    xlabel('time (sec)');
    ylabel(['BOLD signal (' fmriUnits ')']);
%     drawPublishAxis('xLabelOffset', -9/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 4/64);

end


%% FFT of average run, for s0008
% ymax=0;
for r = 1:length(ROIs)
    iROI = ROIs(r);
%     subplot(rows,cols,(r-1)*cols+2)
    subplot(rows,cols,subplots{2})
    for rwd=1:2
        subMeanRun(iRoi,rwd,:) = squeeze(subMeanRunTC(iSub,iRoi,rwd,:));
        subMeanRun(iRoi,rwd,:) = subMeanRun(iRoi,rwd,:) - mean(subMeanRun(iRoi,rwd,:));
        %         subMeanRun(iRoi,rwd,:) = zscore(subMeanRun(iRoi,rwd,:));
        runFFT = squeeze(abs(fft(subMeanRun(iRoi,rwd,:))));
 
        nframes = length(runFFT);
        %from plotMeanFourierAmp.m
        runFFT = runFFT / (nframes/2);
        frequencies = [0:nframes-1]/(nframes*TR);
        runFFT = runFFT(1:floor(nframes/2));
        frequencies = frequencies(1:floor(nframes/2));
        
        plot(frequencies, runFFT,'.-','Color',plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);

        hold on
    end
   xlabel('frequency (cycles/run)');
   ylabel(['Fourier amplitude  (' fmriUnits ')']);

end

%% bar plot of subjects' amplitudes
lineLength = 0.2;
lineWidth = 2;
markerSize = 10;

subplot(rows,cols,subplots{3})
for r = 1:length(ROIs)
    iROI = ROIs(r);
    clear smallSubAmp
    for i=1:length(goodSubs)
        iSub = goodSubs(i);
        for rwd=1:2
            smallSubAmp(i,rwd) = std(squeeze(subResponse(iSub,iRoi,rwd,:)));
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
    
    
    l = size(scatterCmap,1);
    for i=1:length(goodSubs)
%         subColor = jetMap(ceil(i*l/length(goodSubs)),:);
        subColor(i,:) = scatterCmap(1+floor((smallSubDiff(i) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
        plot(1:rewards,smallSubAmp(i,:),'Color',subColor(i,:),'linewidth',linewidth);
        hold on
    end
%     scatter(rwdNum(:),smallSubAmp(:),markerSize, subNum(:),'filled');
    scatter(rwdNum(:),smallSubAmp(:),markerSize, [subColor; subColor] ,'filled');
    
    % MEAN
    for i=1:rewards
%         line([i-lineLength/2 i+lineLength/2], [mean(smallSubAmp(:,i)) mean(smallSubAmp(:,i))], 'color','k','linewidth',lineWidth);
        scatter(i,mean(smallSubAmp(:,i)),markerSize*4,[0 0 0],'filled','d');
        plot(1:rewards,mean(smallSubAmp),'Color',[0 0 0],'linewidth', 2*linewidth);
    end
    xlabel('reward');
    ylabel('response amplitude (std)');
 
end

%%
%% Bar plot of all subjects FFT of average run
% ymax=0;
clear smallSubAmp
subplot(rows,cols,subplots{4})
for r = 1:length(ROIs)
    iROI = ROIs(r);
    for i=1:length(goodSubs)
        iSub = goodSubs(i);
        for rwd=1:2
            subMeanRun(iRoi,rwd,:) = squeeze(subMeanRunTC(iSub,iRoi,rwd,:));
            subMeanRun(iRoi,rwd,:) = subMeanRun(iRoi,rwd,:) - mean(subMeanRun(iRoi,rwd,:));
            %         subMeanRun(iRoi,rwd,:) = zscore(subMeanRun(iRoi,rwd,:));
            runFFT = squeeze(abs(fft(subMeanRun(iRoi,rwd,:))));
            
            nframes = length(runFFT);
            %from plotMeanFourierAmp.m
            runFFT = runFFT / (nframes/2);
            frequencies = [0:nframes-1]/(nframes*TR);
            runFFT = runFFT(1:floor(nframes/2));
            frequencies = frequencies(1:floor(nframes/2));
            rwdFftAmp(i, rwd) = sum(runFFT(ntrials+1:ntrials:end));
            
        end
    end
    fftSubDiff = rwdFftAmp(:, 1) - rwdFftAmp(:, 2);
    minSubFftDiff = min(fftSubDiff);
    maxSubFftDiff = max(fftSubDiff);
    for i=1:length(goodSubs)
%         subColor(i,:) = scatterCmap(1+floor((fftSubDiff(i) - minSubFftDiff)*(l-1)/(maxSubFftDiff-minSubFftDiff)),:);
%keeping the same color from the amplitude difference, previous subplot
        plot(1:rewards,rwdFftAmp(i,:),'Color',subColor(i,:),'linewidth',linewidth);
        hold on
    end
    scatter(rwdNum(:),rwdFftAmp(:),markerSize, [subColor; subColor] ,'filled');

    % MEAN
    for i=1:rewards
        scatter(i,mean(rwdFftAmp(:,i)),markerSize*4,[0 0 0],'filled','d');
        plot(1:rewards,mean(rwdFftAmp),'Color',[0 0 0],'linewidth', 2*linewidth);
    end
    xlabel('reward');
%     ylabel('response amplitude (std)');
    ylabel(['Fourier amplitude  (' fmriUnits ')']);

end


%%
for isubplot=1:length(subplots)
   subplot(rows,cols,subplots{isubplot})

      if isubplot<3
          drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',0,...
              'labelFontSize',7);
          axis square
      else
          drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -50/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
          'xTick',[1 2], 'xTickLabel', {'high','low'},...
            'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',7,...
            'yAxisMajorTickLen',-4/32);
      end
end
set(gcf,'position',[10 10 25 6]);
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig3_' ConcatProjStr '.pdf']);
%%
toc



