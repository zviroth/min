
clear dataFolder s rwdPupil meanPupil stdRwd numRuns numTrials trialLength runSize subPupil
dataFolder{1} = '/Volumes/MH02086153MACDT-Drobo/allMinBehavioral/';
saveFolder = '/Volumes/MH02086153MACDT-Drobo/allMinBehavioral/';

load([saveFolder 'behavioralData.mat'],'dataFolder', 'subdirs', 'subPupil', 'runSize', 'rwdPupil','meanPupil','stdRwd','diffPupil','expName');

%%

figure(1)
clf
rows = 2;
cols = ceil(length(subdirs)/rows);
for iSub = 1:length(subdirs)
    subplot(rows,cols,iSub)
    for rwd=1:2%H=1, L=2
        plot(meanPupil{iSub,rwd})
        hold all
    end
    title(expName{iSub,1})
end


figure(2)
clf
for iSub = 1:length(subdirs)
    minlength(iSub) = min(length(meanPupil{iSub,1}), length(meanPupil{iSub,2}));
    diffPupil{iSub} = meanPupil{iSub,1}(1:minlength(iSub)) - meanPupil{iSub,2}(1:minlength(iSub));
    subplot(rows,cols,iSub)
    plot(diffPupil{iSub});
    hold on
    plot(zeros(1,length(diffPupil{iSub})),'k')
    title(expName{iSub,1})
     
end

figure(3)
clf
rows=1;
cols=2;
subplot(rows,cols,1)
allMin = min(minlength);
clear allSubPupil
for iSub = 1:length(subdirs)
    for rwd=1:2%H=1, L=2
        allSubPupil(iSub,rwd,:) = meanPupil{iSub,rwd}(1:allMin);
    end
end
groupPupil = squeeze(mean(allSubPupil));
plot(groupPupil','linewidth',2);

% for rwd=1:2
%    plot(groupPupil(rwd,:));
%    hold all
% end
subplot(rows,cols,2)
plot(groupPupil(1,:)-groupPupil(2,:),'k','linewidth',2)