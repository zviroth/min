dataFolder = '/Volumes/MH02086153MACDT-Drobo/decodingAnalysis/rwd/';
onlyCorrect=0;
onlyCorrectString = '';
if onlyCorrect
    onlyCorrectString = '_correct';
end
load([dataFolder 'subMeanRoiTC' onlyCorrectString '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', 'meanResponse', 'stdResponse',...
    'roiTC', 'allTrials', ...
    'subdirs', 'roiNames','subTrialResponse','subRunResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect');
roiNames(5:end) = [];
roiNames([1:2 4:end]) = [];
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
clear subMeanResponse trialStd meanTrialStd runStd meanRunStd subRwdStd trialFFTamp trialFFTphase meanTrialFFTamp meanTrialFFTphase 
clear stdTrialFFTphase runFFTamp runFFTphase meanRunFFTamp meanRunFFTphase stdRunFFTphase subFFTamp subFFTphase
clear trialMinMax meanTrialMinMax runMinMax meanRunMinMax subRwdMinMax
for iSub = 1:length(subdirs)
    subplot(rows,cols,iSub)
    for iRoi= 1:length(roiNames)
        subMeanResponse{iSub,iRoi} = mean(squeeze(subResponse(iSub,iRoi,:,:)));%averaged over reward type
        for rwd=1:2
            
%             STD
            trialStd{iSub,iRoi,rwd} = std(subTrialResponse{iSub,iRoi,rwd});%std per trial
            meanTrialStd(iSub,iRoi,rwd) = mean(trialStd{iSub,iRoi,rwd});
            runStd{iSub,iRoi,rwd} = std(subRunResponse{iSub,iRoi,rwd});%std per run
            meanRunStd(iSub,iRoi,rwd) = mean(runStd{iSub,iRoi,rwd});
            subRwdStd(iSub,iRoi,rwd) = std(squeeze(subResponse(iSub,iRoi,rwd,:)));
            
%             FFT
            trialFFT = fft(subTrialResponse{iSub,iRoi,rwd});
            trialFFTamp{iSub,iRoi,rwd} = abs(trialFFT(2,:));
            trialFFTphase{iSub,iRoi,rwd} = angle(trialFFT(2,:));
            meanTrialFFTamp(iSub,iRoi,rwd) = mean(trialFFTamp{iSub,iRoi,rwd});
%             meanTrialFFTphase(iSub,iRoi,rwd) = mean(trialFFTphase{iSub,iRoi,rwd});
%             stdTrialFFTphase(iSub,iRoi,rwd) = std(trialFFTphase{iSub,iRoi,rwd});
            [s s0] =  circ_std(trialFFTphase{iSub,iRoi,rwd}');
            stdTrialFFTphase(iSub,iRoi,rwd) = s0;
            
            runFFT = fft(subRunResponse{iSub,iRoi,rwd});
            runFFTamp{iSub,iRoi,rwd} = abs(runFFT(2,:));
            runFFTphase{iSub,iRoi,rwd} = angle(runFFT(2,:));
            meanRunFFTamp(iSub,iRoi,rwd) = mean(runFFTamp{iSub,iRoi,rwd});
%             meanRunFFTphase(iSub,iRoi,rwd) = mean(runFFTphase{iSub,iRoi,rwd});
%             stdRunFFTphase(iSub,iRoi,rwd) = std(runFFTphase{iSub,iRoi,rwd});
            [s s0] =   circ_std(runFFTphase{iSub,iRoi,rwd}');
            stdRunFFTphase(iSub,iRoi,rwd) = s0;
            
            subFFT = fft(squeeze(subResponse(iSub,iRoi,rwd,:)));
            subFFTamp(iSub,iRoi,rwd) = abs(subFFT(2));
            subFFTphase(iSub,iRoi,rwd) = angle(subFFT(2));
            
            %MIN-MAX
            trialMinMax{iSub,iRoi,rwd} = max(subTrialResponse{iSub,iRoi,rwd})-min(subTrialResponse{iSub,iRoi,rwd});% per trial
            meanTrialMinMax(iSub,iRoi,rwd) = mean(trialMinMax{iSub,iRoi,rwd});
            runMinMax{iSub,iRoi,rwd} = max(subRunResponse{iSub,iRoi,rwd})-min(subRunResponse{iSub,iRoi,rwd});% per run
            meanRunMinMax(iSub,iRoi,rwd) = mean(runMinMax{iSub,iRoi,rwd});
            subRwdMinMax(iSub,iRoi,rwd) = max(squeeze(subResponse(iSub,iRoi,rwd,:))) - min(squeeze(subResponse(iSub,iRoi,rwd,:)));
            
            % plot mean timecourse
            plot(squeeze(subResponse(iSub,iRoi,rwd,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
            title(subdirs(iSub).name(1:5));
            hold on
        end
    end
%     %correlation between V1-DMN amplitudes 
%     for rwd=1:2
%         trialStdCorr{iSub} = corr(trialStd{iSub,1,rwd}', trialStd{iSub,2,rwd}');
%         runStdCorr{iSub} = corr(runStd{iSub,1,rwd}', runStd{iSub,2,rwd}');
%     end
end
legendCellArray = [];
for iRoi= 1:length(roiNames)
    legendCellArray{(iRoi-1)*2+1} = [roiNames{iRoi} ': H'];
    legendCellArray{iRoi*2} = [roiNames{iRoi} ': L'];
end
legend(legendCellArray);
set(gcf,'position',[250 500 1200 400]);

%% DIFFERENT MEASURES OF RESPONSE, PER SUBJECT/ROI
figure(2)
clf
rows=3;
cols=3;
% STD
subplot(rows,cols,1)
plot(meanTrialStd(:,:,1)-meanTrialStd(:,:,2), 'linewidth', 1); title('trials STD');
subplot(rows,cols,2)
plot(meanRunStd(:,:,1)-meanRunStd(:,:,2), 'linewidth', 1); title('runs STD');
subplot(rows,cols,3)
plot(subRwdStd(:,:,1)-subRwdStd(:,:,2), 'linewidth', 1); title('mean STD');
% FFT
phaseNorm = 0.01;
subplot(rows,cols,4)
plot(meanTrialFFTamp(:,:,1)-meanTrialFFTamp(:,:,2), 'linewidth', 1); title('trials FFT');
hold on
plot(phaseNorm*(stdTrialFFTphase(:,1,1)-stdTrialFFTphase(:,1,2)),'b.', 'linewidth', 1);
corr(meanTrialFFTamp(:,1,1)-meanTrialFFTamp(:,1,2), stdTrialFFTphase(:,1,2)-stdTrialFFTphase(:,1,1))
subplot(rows,cols,5)
plot(meanRunFFTamp(:,:,1)-meanRunFFTamp(:,:,2), 'linewidth', 1); title('runs FFT');
hold on
plot(phaseNorm*(stdRunFFTphase(:,1,1)-stdRunFFTphase(:,1,2)),'b.', 'linewidth', 1);
corr(meanRunFFTamp(:,1,1)-meanRunFFTamp(:,1,2),stdRunFFTphase(:,1,2)-stdRunFFTphase(:,1,1))
subplot(rows,cols,6)
plot(subFFTamp(:,:,1)-subFFTamp(:,:,2), 'linewidth', 1); title('mean FFT');
% MIN-MAX
subplot(rows,cols,7)
plot(meanTrialMinMax(:,:,1)-meanTrialMinMax(:,:,2), 'linewidth', 1); title('trials max-min');
subplot(rows,cols,8)
plot(meanRunMinMax(:,:,1)-meanRunMinMax(:,:,2), 'linewidth', 1); title('runs max-min');
subplot(rows,cols,9)
plot(subRwdMinMax(:,:,1)-subRwdMinMax(:,:,2), 'linewidth', 1); title('mean max-min');

% add zero baseline lines
for r=1:rows
    for c=1:cols
        subplot(rows,cols,c + (r-1)*cols)
        hold on
        plot(1:length(subdirs),zeros(1,length(subdirs)),'k--');
%         yticks([]);
        xlim([1 length(subdirs)])
    end
end
set(gcf,'position',[200 150 800 600]);
%% BAR PLOTS
clear meanRwdStd
figure(3)
clf
xshift=0.2;
scatterSize = 30;
linewidth = 3;
linelength = 0.2;
%line per subject per ROI, connecting low and high reward
for iRoi= 1:length(roiNames)
    line([iRoi-0.5*xshift iRoi+0.5*xshift], [subRwdStd(:,iRoi,1) subRwdStd(:,iRoi,2)],'linewidth',linewidth/2);
    hold all
    pval(iRoi) = ranksum(subRwdStd(:,iRoi,1), subRwdStd(:,iRoi,2));
end
%black bars for average across subjects
for rwd=1:2
    meanRwdStd(:,rwd) = mean(subRwdStd(:,:,rwd)); %per ROI
    for iRoi= 1:length(roiNames)
        scatterCenter = iRoi+(rwd-1.5)*xshift;
        scatter(scatterCenter*ones(length(subdirs),1), subRwdStd(:,iRoi,rwd),scatterSize, plotColors{rwd});
        hold all
        line([scatterCenter-linelength/2 scatterCenter+linelength/2], [meanRwdStd(iRoi,rwd) meanRwdStd(iRoi,rwd)], 'color','k','linewidth',linewidth);
    end 
end
%black line connecting between means
for iRoi= 1:length(roiNames)
    line([iRoi-0.5*xshift iRoi+0.5*xshift], [meanRwdStd(iRoi,1) meanRwdStd(iRoi,2)],'linewidth',linewidth/2, 'color', 'k','linewidth',linewidth);
end


xlim([0 length(roiNames)+1]);
xticks(1:length(roiNames));
xticklabels(roiNames);


% %%
% figure(4)
% clf
% cols=ceil(length(subdirs)/2);
% rows = ceil(length(subdirs)/cols);
% for iSub = 1:length(subdirs)
%     subplot(rows,cols,iSub)
%     for iRoi= 1:2%length(roiNames)
%         for rwd=1:2
%             plot(squeeze(subResponse(iSub,iRoi,rwd,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
%             hold on
%         end
%     end
% end

%%
figure(5)
for iRoi= 1:length(roiNames)
    for rwd=1:2
        plot(squeeze(meanResponse(rwd,iRoi,:)),plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
        hold on
    end
end
legend(legendCellArray);
%%
figure(6)
clf
for iRoi= 1:length(roiNames)
    for rwd=1:2
        stdAmplitude(iRoi,rwd) = std(mean(allTrials{iRoi,rwd},2));
        plot(mean(allTrials{iRoi,rwd},2),plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
%         errorbar(1:size(allTrials{iRoi,rwd},1),std(allTrials{iRoi,rwd},1,2), 'Color', plotColors{rwd});
        errorbar(mean(allTrials{iRoi,rwd},2),std(allTrials{iRoi,rwd},1,2), 'Color', plotColors{rwd});
%         plot(squeeze(meanResponse(rwd,iRoi,:)),plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
        hold on
    end
end


% %PERMUTATIONS
% nperms=100000;
% for rwd=1:2
%     numTrials(rwd) = size(allTrials{1,rwd},2);
%     if rwd==1
%         firstTrial(rwd)=1;
%     else %rwd==2
%         firstTrial(rwd)=numTrials(1)+1;
%     end
% end
% for iRoi=1:length(roiNames)
%    roiTrials{iRoi} = [allTrials{iRoi,1} allTrials{iRoi,2}]; 
% end
% 
% for p=1:nperms
%     randOrder = randperm(numTrials(1)+numTrials(2));
%     for iRoi= 1:length(roiNames)
%         for rwd=1:2
%             permTrials = roiTrials{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(rwd)-1));
%             permStdAmplitude(iRoi,rwd,p) = std(mean(permTrials,2));
%         end
%     end
% end
% for iRoi= 1:length(roiNames)
%     realDiff = stdAmplitude(iRoi,1) - stdAmplitude(iRoi,2);
%     permDiff = squeeze(permStdAmplitude(iRoi,1,:) - permStdAmplitude(iRoi,2,:));
%     pVal(iRoi) = sum(permDiff>realDiff)/nperms;
% end
%     
% pVal
%%
for iSub = 1:length(subdirs)
   for rwd=1:2
       meanPropCorrect(iSub,rwd) = mean(propCorrect{iSub,rwd});
   end
end
