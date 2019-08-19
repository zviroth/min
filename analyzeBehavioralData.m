saveFolder = '/Volumes/MH02086153MACDT-Drobo/allMinBehavioral/';

load([saveFolder 'behavioralData.mat'],'dataFolder', 'subFolders', 'subPupil', 'runSize', 'rwdPupil',...
    'meanPupil','stdRwd','diffPupil','expName',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh','rwdLevel','numTrials');

plotColors = { [0 0 1],[1 0 0],[0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
%%
i=0;

i=i+1;
figure(i)
clf
rows = 2;
cols = ceil(length(subFolders)/rows);
minLength=Inf;
for iSub = 1:length(subFolders)
    subplot(rows,cols,iSub)
    for rwd=1:2
        plot(meanPupil{iSub,rwd})
        hold all
        if minLength>length(meanPupil{iSub,rwd})
            minLength = length(meanPupil{iSub,rwd});
        end
    end
    title(expName{iSub,1})
end
minLength = 2000;
for rwd=1:2
    groupAll{rwd}(1,:) = meanPupil{1,rwd}(1:minLength);
    for iSub=2:length(subFolders)
        groupAll{rwd}(iSub,:) =  meanPupil{iSub,rwd}(1:minLength);
    end
    groupMean{rwd} = mean(groupAll{rwd});
    groupStd{rwd} = std(groupAll{rwd});
end

i=i+1;
figure(i)
clf
clear diffPupil
for iSub = 1:length(subFolders)
%     minlength = min(length(meanPupil{iSub,1}), length(meanPupil{iSub,2}));
    diffPupil(iSub,:) = meanPupil{iSub,1}(1:minLength) - meanPupil{iSub,2}(1:minLength);
    subplot(rows,cols,iSub)
    plot(diffPupil(iSub,:));
    hold on
    plot(zeros(1,length(diffPupil(iSub,:))),'k')
    title(expName{iSub,1})
end


dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;
%%
i=i+1;
figure(i)
clf
subplot(1,2,1)
for rwd=1:2
    dsErrorsurface((1:minLength)*2, groupMean{rwd}(1:minLength), groupStd{rwd}(1:minLength)./sqrt(length(subFolders)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
end
for rwd=1:2
    plot((1:minLength)*2,groupMean{rwd}(1:minLength),'Color',plotColors{rwd},'linewidth',linewidth);
end
title('group mean');
axis square
box on
xlabel('time (ms)');
ylabel('pupil size (arb. units)');
% legend({'high','low'},'location','northeast');
drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -10/64, 'xAxisMargin', 4/64, 'yAxisMargin', 4/64);
% legend('H','L','location','northeast');
subplot(1,2,2)
meanDiff = mean(diffPupil);
stdDiff = std(diffPupil);
dsErrorsurface((1:minLength)*2, meanDiff(1:minLength), stdDiff(1:minLength)./sqrt(length(subFolders)), 'k',dsSurfaceAlpha);
hold on
plot((1:minLength)*2, meanDiff(1:minLength),'k','linewidth',linewidth);
title(' H - L');
axis square
box on
xlabel('time (ms)');
ylabel('\Delta pupil size (arb. units)');
drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -8/64, 'xAxisMargin', 4/64, 'yAxisMargin', 2/64);

%%
i=i+1;
figure(i)
clf
iSub=2;
subplot(1,2,1)
for rwd=1:2
    dsErrorsurface((1:minLength)*2, meanPupil{iSub,rwd}(1:minLength), stdRwd{iSub,rwd}(1:minLength)./sqrt(numTrials(iSub,rwd)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
end
for rwd=1:2
    plot((1:minLength)*2,meanPupil{iSub,rwd}(1:minLength),'Color',plotColors{rwd},'linewidth',linewidth);
end
title(subFolders{iSub}(1:5));
axis square
box on
xlabel('time (ms)');
ylabel('pupil size (arb. units)');
% legend({'high','low'},'location','northeast');
drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -10/64, 'xAxisMargin', 4/64, 'yAxisMargin', 4/64);
% legend('H','L','location','northeast');

subplot(1,2,2)
plot((1:minLength)*2, diffPupil(iSub,:),'k','linewidth',linewidth);
title([subFolders{iSub}(1:5) ' H - L']);
axis square
box on
xlabel('time (ms)');
ylabel('\Delta pupil size (arb. units)');
drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -8/64, 'xAxisMargin', 4/64, 'yAxisMargin', 4/64);
%%
mean(subMeanCorrectness)
mean(subMeanRT)
mean(subMedianRT)
mean(subMeanThresh)