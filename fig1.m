close all;
clear all;
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
plotColors = { [0 0 1],[1 0 0],[0 1 0], [0.5 1 0.2]};
plotColors = { [1 0 0],[0 0 1],[0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
fontsize=9;
dsSurfaceContrast = 1;
dsSurfaceAlpha = 0.15;
%%
i=0;
minLength = 2000;
for rwd=1:2
    groupAll{rwd}(1,:) = meanPupil{1,rwd}(1:minLength);
    for iSub=2:length(subFolders)
        groupAll{rwd}(iSub,:) =  meanPupil{iSub,rwd}(1:minLength);
    end
    groupMean{rwd} = mean(groupAll{rwd});
    groupStd{rwd} = std(groupAll{rwd});
end
clear diffPupil groupAll
for iSub = 1:length(subFolders)
%     minlength = min(length(meanPupil{iSub,1}), length(meanPupil{iSub,2}));
    diffPupil(iSub,:) = meanPupil{iSub,1}(1:minLength) - meanPupil{iSub,2}(1:minLength);
end


% i=i+1;
% figure(i)
% clf
% rows=2;
% cols=ceil(length(subFolders)/2);
% for iSub=1:length(subFolders)
%     subplot(rows,cols,iSub)
%     for rwd=1:2
%         dsErrorsurface((1:minLength)*2, meanPupil{iSub,rwd}(1:minLength), stdRwd{iSub,rwd}(1:minLength)./sqrt(numTrials(iSub,rwd)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
%         hold on
%     end
%     for rwd=1:2
%         plot((1:minLength)*2,meanPupil{iSub,rwd}(1:minLength),'Color',plotColors{rwd},'linewidth',linewidth);
%     end
% end

% i=i+1;
% figure(i)
% clf
% for iSub=1:length(subFolders)
%     subplot(rows,cols,iSub)
%     
%     plot((1:minLength)*2, diffPupil(iSub,:),'k','linewidth',linewidth);
% end


%%
i=i+1;
figure(i)
clf
rows=2;
cols=2;

%%

iSub=6;
subplot(rows,cols,1)
for rwd=1:2
    dsErrorsurface((1:minLength)*2, meanPupil{iSub,rwd}(1:minLength), stdRwd{iSub,rwd}(1:minLength)./sqrt(numTrials(iSub,rwd)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
end
for rwd=1:2
    plot((1:minLength)*2,meanPupil{iSub,rwd}(1:minLength),'Color',plotColors{rwd},'linewidth',linewidth);
end

subplot(rows,cols,2)
plot((1:minLength)*2, diffPupil(iSub,:),'k','linewidth',linewidth);
hline(0);

subplot(rows,cols,3)
for rwd=1:2
    dsErrorsurface((1:minLength)*2, groupMean{rwd}(1:minLength), groupStd{rwd}(1:minLength)./sqrt(length(subFolders)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
end
for rwd=1:2
    plot((1:minLength)*2,groupMean{rwd}(1:minLength),'Color',plotColors{rwd},'linewidth',linewidth);
end

subplot(rows,cols,4)
meanDiff = mean(diffPupil);
stdDiff = std(diffPupil);
dsErrorsurface((1:minLength)*2, meanDiff(1:minLength), stdDiff(1:minLength)./sqrt(length(subFolders)), 'k',dsSurfaceAlpha);
hold on
plot((1:minLength)*2, meanDiff(1:minLength),'k','linewidth',linewidth);
hline(0);

for isubplot=1:rows*cols
   subplot(rows,cols,isubplot)
   xlabel('time (ms)');
   axis square
   if mod(isubplot,2)>0
       ylabel('pupil size (arb. units)');
        drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 4/64,'labelFontSize',fontsize);
   else
       ylabel('\Delta pupil size (arb. units)');
       drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 2/64,'labelFontSize',fontsize);
   end
end
set(gcf,'position',[10 10 18 15]);

print('-painters','-dpdf','~/Documents/MATLAB/min/figures/fig1_pupil.pdf');

%%
mean(subMeanCorrectness)
mean(subMeanRT)
mean(subMedianRT)
mean(subMeanThresh)

%% mean + STD of thresholds
subThresh = mean(subMeanThresh,2);%mean over rwd, which is already averaged over staircases and runs
mean(subThresh)
std(subThresh)
% keyboard
%%

% %% A figure just for the legend!
% % plotColors = { [0.9 0.1 0.1],[0.1 0.1 0.9],[0 1 0], [0.5 1 0.2]};
% i=i+1;
% figure(i)
% clf
% for rwd=1:2
%    plot(rwd+[1:10], 'Color', plotColors{rwd},'linewidth',2);
%    hold on
% end
% [~, hobj, ~, ~] = legend('high reward', 'low reward');
% hl = findobj(hobj,'type','line');
% set(hl,'LineWidth',1.5);
% ht = findobj(hobj,'type','text');
% set(ht,'FontSize',12);
% print('-painters','-dpdf','~/Documents/MATLAB/min/figures/legend.pdf');
%%
% print('-painters','-dpdf','~/Documents/MATLAB/min/figures/fig1_pupil.pdf');