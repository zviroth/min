close all
clear all
tic
saveFolder = '/Volumes/MH02086153MACDT-Drobo/allMinBehavioral/';
dataFolder = saveFolder;
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};

load([dataFolder 'eyeData.mat'], 'dataFolder', 'subFolders', 'samplerate',  ...
    'numSubs', 'onlyCorrect', ...
    'sacRwd','binnedSac', 'smoothSac', ...
    'rwdPupil','meanPupil',...
    'sacRun','pupilRun', 'numRuns','startTimes','endTimes','eyetrackerTime',...
    'listRwdRuns', 'concatRwdTypeNum','expName',...
    'sacRunSmooth', 'sacRunSmoothTrialZeroFilled','sacRunSmoothTrial','L');

clear sacPupilFunc sacPupilFuncRun resid residRun subMeanPupilLength sacMat sacPupilRun
% numSubs=2;
% numSubs=3;
% sacRun{iSub,rwdTypeNum,r}; %saccade onset
% pupilRun{iSub,rwdTypeNum,r}; %pupil size
deconvLength= 500; %in eyetracker samples
for iSub=1:numSubs
%     numRuns(iSub) = size(pupilRun,3);%what if there is a different number of high and low runs?
    concatSac=[];
    concatPupil=[];
    concatSacMat=[];
    for rwd=1:2

        for r=1:numRuns(iSub,rwd)

%             %detrend, convert to % signal change
%             goodtimepoints = ~isnan(pupilRun{iSub,rwd,r});
%             temp = pupilRun{iSub,rwd,r}(goodtimepoints);
%             temp = percentTSeries(temp')';
%             pupilRun{iSub,rwd,r}(goodtimepoints) = temp;


            runLength = length(pupilRun{iSub,rwd,r}(:));
            sacVec = zeros(1,runLength);
            sacVec(sacRun{iSub,rwd,r}) = ones;%one for each microsaccade onset
            
            sacMat{iSub,rwd,r} = zeros(deconvLength,runLength);
            for j=1:deconvLength
                sacMat{iSub,rwd,r}(j,:) = [zeros(1,j-1) sacVec(1:end-j+1)];%shift and pad with zeros
            end
            sacMat{iSub,rwd,r}(deconvLength+1,:) = ones;%constant predictor
            
            goodTimepoints = ~isnan(pupilRun{iSub,rwd,r});
            %single run pupil function, for illustration purposes...
            sacPupilFuncRun{iSub}(rwd,r,:) = sacMat{iSub,rwd,r}(:,goodTimepoints)'\pupilRun{iSub,rwd,r}(goodTimepoints)';%deconvolution
%             residRun{iSub,rwd,r}(:) = pupilRun{iSub,rwd,r}' - sacMat'*squeeze(sacPupilFuncRun{iSub}(rwd,r,:));%residual timecourse
            
            concatSacMat = [concatSacMat sacMat{iSub,rwd,r}];
            concatSac = [concatSac; sacRun{iSub,rwd,r}+length(concatPupil)];
            concatPupil = [concatPupil pupilRun{iSub,rwd,r}];
            
        end
    end
    goodTimepoints = ~isnan(concatPupil);
    sacPupilFunc(iSub,:) = concatSacMat(:,goodTimepoints)'\concatPupil(goodTimepoints)';
    
    %set the function integral to zero
    offset = sum(sacPupilFunc(iSub,1:end-1))./length(sacPupilFunc(iSub,1:end-1));
    sacPupilFunc(iSub,1:end-1) = sacPupilFunc(iSub,1:end-1)-offset;
    
%      sacPupilFunc(iSub,end) =  sacPupilFunc(iSub,end) + offset;
    
%     sacPupilFunc(iSub,:) = sacPupilFunc(iSub,:) - sacPupilFunc(iSub,1);%subtract the first time point of the pupil-ms function from the entire function, so that it starts at zero
    
    
%     resid{iSub} = concatPupil' - concatSacMat(1:end-1,:)'*sacPupilFunc(iSub,1:end-1)';
%     resid{iSub} = resid{iSub} + sacPupilFunc(iSub,end);%add the constant
    
    %perform regression for each run, using the concatenated data's ms
    %pupil function
    for rwd=1:2
        for r=1:numRuns(iSub,rwd)
            sacPupilRun{iSub,rwd,r} = sacMat{iSub,rwd,r}(1:end-1,:)'*squeeze(sacPupilFunc(iSub,1:end-1))';%ms pupil function convolved with all microsaccade onsets
            residRun{iSub,rwd,r}(:) = pupilRun{iSub,rwd,r}' - sacMat{iSub,rwd,r}(1:end-1,:)'*squeeze(sacPupilFunc(iSub,1:end-1))';%residual timecourse. 
%             residRun{iSub,rwd,r}(:) = residRun{iSub,rwd,r}(:) + sacPupilFunc(iSub,end); %add the constant 
        end
    end

    %cut residuals into trials, and average
    for rwd=1:2
        residRwdPupil{iSub,rwd} = NaN(size(rwdPupil{iSub,rwd}));%numTrials X maxTrialLength
        regressRwdPupil{iSub,rwd} = NaN(size(rwdPupil{iSub,rwd}));
        trialCounter=0;%trial counter across runs
       %get residual for each run
       for r=1:numRuns(iSub,rwd)
           for itrial=1:size(startTimes{iSub,rwd,r},2)
               
               thisTrialTimes = startTimes{iSub,rwd,r}(itrial):endTimes{iSub,rwd,r}(itrial);
               [~, t1] = min(abs(eyetrackerTime{iSub,rwd,r}(1:length(pupilRun{iSub,rwd,r}))-thisTrialTimes(1)));
               t1=t1(1);
               [~, t2] = min(abs(eyetrackerTime{iSub,rwd,r}(1:length(pupilRun{iSub,rwd,r}))-thisTrialTimes(end)));
               t2=t2(1);
               
               trialResidual = residRun{iSub,rwd,r}(t1:t2);
%                trialResidual = pupilRun{iSub,rwd,r}(t1:t2); % TESTING

               
               
%                trialPupil = residRun{iSub,rwd,r}(startTimes{iSub,rwd,r}(itrial):endTimes{iSub,rwd,r}(itrial)-1);%endTime is startTime of next trial
               trialCounter = trialCounter+1;
               residRwdPupil{iSub,rwd}(trialCounter,1:length(trialResidual)) = trialResidual;
               
               trialRegressed = sacPupilRun{iSub,rwd,r}(t1:t2);
               regressRwdPupil{iSub,rwd}(trialCounter,1:length(trialResidual)) = trialRegressed;%this is regressed out, the pupil modulation from microsaccades
           end
       end
       
       %average across trials
       residMeanPupil{iSub,rwd} = nanmean(residRwdPupil{iSub,rwd});
       regressMeanPupil{iSub,rwd} = nanmean(regressRwdPupil{iSub,rwd});
    end
    
    

    l = min(length(meanPupil{iSub,1}), length(meanPupil{iSub,2}));
    subMeanPupilLength(iSub) = l;
end

%%
allMeanPupil = zeros(2,min(subMeanPupilLength));
allResidMeanPupil = zeros(2,min(subMeanPupilLength));
sumBinnedSac = zeros(2,min(subMeanPupilLength));
sumSmoothSac = zeros(2,min(subMeanPupilLength));
allMeanStimPupil = zeros(2,min(subMeanPupilLength));
allResidStimPupil = zeros(2,min(subMeanPupilLength));
allMeanNullPupil = zeros(2,min(subMeanPupilLength));
allResidNullPupil = zeros(2,min(subMeanPupilLength));
allRegressMeanPupil = zeros(2,min(subMeanPupilLength));

for iSub=1:numSubs
    for rwd=1:2
        allMeanPupil(rwd,:) = allMeanPupil(rwd,:) + meanPupil{iSub,rwd}(1:min(subMeanPupilLength))';
        allResidMeanPupil(rwd,:) = allResidMeanPupil(rwd,:) + residMeanPupil{iSub,rwd}(1:min(subMeanPupilLength));
        
        sumBinnedSac(rwd,:) =  sumBinnedSac(rwd,:) + binnedSac{iSub,rwd}(1:min(subMeanPupilLength));
        sumSmoothSac(rwd,:) =  sumSmoothSac(rwd,:) + smoothSac{iSub,rwd}(1:min(subMeanPupilLength));
          
        allRegressMeanPupil(rwd,:) = allRegressMeanPupil(rwd,:) + regressMeanPupil{iSub,rwd}(1:min(subMeanPupilLength));
    end
end
allMeanPupil = allMeanPupil./numSubs;
allResidMeanPupil = allResidMeanPupil./numSubs;
sumBinnedSac = sumBinnedSac./numSubs;
sumSmoothSac = sumSmoothSac./numSubs;
        
allRegressMeanPupil = allRegressMeanPupil./numSubs;

%%
rwdTypes = {'high','low'};
figure(1); clf%pupil size

rows=4;
cols=numSubs;
for iSub=1:numSubs
    minPup = min(min(meanPupil{iSub,1}),min(meanPupil{iSub,2})) - 100;
    maxPup = max(max(meanPupil{iSub,1}),max(meanPupil{iSub,2})) + 100;
    for rwd=1:2
        subplot(rows,cols,iSub + (rwd-1)*cols)
        plot(rwdPupil{iSub,rwd}')
        hold all
        plot(meanPupil{iSub,rwd}, 'color', plotColors{rwd}, 'linewidth', 2)
        
        titleString = [rwdTypes{rwd} ' reward: ' num2str(size(rwdPupil{iSub,rwd},1)) ' trials'];
        title(titleString);
%             ylim([minPup maxPup]);
        xlabel('time (ms)')
        ylabel('pupil size (pixels)');
            xticks(0:1500:3000);
        xticklabels(2*[0:1500:3000]);
        xlim([0 1500]);
%         xticks(0:1000:8000);
%         xticklabels(2*[0:1000:8000]);
    end

    subplot(rows,cols,iSub+2*cols)
    for rwd=1:2
        plot(meanPupil{iSub,rwd}, 'Color', plotColors{rwd}, 'linewidth', 1)
        hold on
    end
    if iSub==1
        legend('high', 'low')
    end
    xlabel('time (ms)')
    ylabel('pupil size (pixels)');
    ylim([minPup maxPup]);
    xticks(0:1000:8000);
    xticklabels(2*[0:1000:8000]);
    subplot(rows,cols,iSub+3*cols)

    l=subMeanPupilLength(iSub)-100;
    diffPupil = meanPupil{iSub,1}(1:l) - meanPupil{iSub,2}(1:l);
    plot(zeros(l,1),'--');
    hold all
    plot(diffPupil, 'color', 'k', 'linewidth', 1)
    title('high minus low');
    xlabel('time (ms)')
    xticks(0:1000:8000);
    xticklabels(2*[0:1000:8000]);
    ylabel('difference in pupil size (pixels)');
    set(gcf, 'position', [130 130 600 1000]);
end

figure(2); clf%microsaccades
rows=2;
cols=numSubs;
% L = 100;
% filter = ones(L,1);
% % filter = [1:L/2 L/2:-1:1];
% filter = filter./sum(filter);
for iSub=1:numSubs
    subplot(rows,cols,iSub)
    for rwd=1:2
        plot(binnedSac{iSub,rwd})
        hold on
    end
    ylim([0 20]);
    xlabel('time (ms)')
    ylabel('# of microsaccade onsets');
    title('microsaccade onsets: original data');
%     xticks(0:4000:8000);
%     xticklabels(2*[0:4000:8000]);
    xticks(0:1500:3000);
    xticklabels(2*[0:1500:3000]);
    xlim([0 1500]);
    subplot(rows,cols,iSub + cols)
    for rwd=1:2
%         smoothSac{iSub,rwdTypeNum} = conv(binnedSac{iSub,rwdTypeNum}, filter, 'same');
        plot(smoothSac{iSub,rwd},'linewidth', 1)
        hold on
    end
        if iSub==1
        legend('high', 'low')
    end
    xlabel('time (ms)')
    ylabel('# of microsaccade onsets');
    title('microsaccade onsets: smoothed data');
%     xticks(0:4000:8000);
%     xticklabels(2*[0:4000:8000]);
        xticks(0:1500:3000);
    xticklabels(2*[0:1500:3000]);
    xlim([0 1500]);
    set(gcf, 'position', [100 100 600 500]);
    
end


%%
figure(3); clf
rows=2;
cols=numSubs;
for iSub=1:numSubs
    subplot(rows,cols,iSub)
    for rwd=1:2
        for r=1:numRuns(iSub)
            plot(squeeze(sacPupilFuncRun{iSub}(rwd,r,1:end-1)) , 'color',plotColors{rwd},'linewidth',0.5);
            hold all
        end
        plot(squeeze(mean(sacPupilFuncRun{iSub}(rwd,:,1:end-1))) , 'color',plotColors{rwd},'linewidth',2);
    end
    temp=squeeze(mean(sacPupilFuncRun{iSub}(:,:,1:end-1)));
    plot(squeeze(mean(temp)) , 'color','k','linewidth',3);
    yl = ylim;
    
    subplot(rows,cols,iSub+cols)
    plot(sacPupilFunc(iSub,1:end-1) , 'linewidth',2);% ignore constant
    ylim(yl);
end

%%
figure(4); clf
rows=2;
cols=2;
subplot(rows,cols,1)
plot(sumBinnedSac');
title('mean microsaccade onsets');
xlim([0 1800]);
subplot(rows,cols,2)
plot(sumSmoothSac');
title('smoothed mean ms rate');
xlim([0 1800]);

subplot(rows,cols,3)
plot(squeeze(nanmean(sacPupilFunc(:,1:end-1))),'linewidth', 3);%ignore constant
title('pupil-microsaccade function');

clear meanConv
for rwd=1:2
    temp = conv(sumSmoothSac(rwd,:)',squeeze(nanmean(sacPupilFunc(:,1:end-1))));
    meanConv(rwd,:) = temp(1:length(sumSmoothSac));
end

subplot(rows,cols,4)
plot(meanConv')
title('pupil-ms function convolved with mean ms rate');
xlim([0 1800]);

%%
rwdTypes = {'high','low'};
figure(5); clf%pupil size

rows=4;
cols=numSubs;
for iSub=1:numSubs
    minPup = min(min(residMeanPupil{iSub,1}),min(residMeanPupil{iSub,2})) - 100;
    maxPup = max(max(residMeanPupil{iSub,1}),max(residMeanPupil{iSub,2})) + 100;
    for rwd=1:2
        subplot(rows,cols,iSub + (rwd-1)*cols)
        plot(residRwdPupil{iSub,rwd}')
        hold all
        plot(residMeanPupil{iSub,rwd}, 'color', plotColors{rwd}, 'linewidth', 2)
        
        titleString = [rwdTypes{rwd} ' reward: ' num2str(size(residRwdPupil{iSub,rwd},1)) ' trials'];
        title(titleString);
%             ylim([minPup maxPup]);
        xlabel('time (ms)')
        ylabel('pupil size (pixels)');
        xticks(0:1000:8000);
        xticklabels(2*[0:1000:8000]);
    end

    subplot(rows,cols,iSub+2*cols)
    for rwd=1:2
        plot(residMeanPupil{iSub,rwd}, 'Color', plotColors{rwd}, 'linewidth', 1)
        hold on
    end
        if iSub==1
        legend('high', 'low')
    end
    xlabel('time (ms)')
    ylabel('pupil size (pixels)');
    ylim([minPup maxPup]);
    xticks(0:1000:8000);
    xticklabels(2*[0:1000:8000]);
    subplot(rows,cols,iSub+3*cols)
    l = min(length(residMeanPupil{iSub,1}), length(residMeanPupil{iSub,2}));
    l=l-100;
    diffPupil = residMeanPupil{iSub,1}(1:l) - residMeanPupil{iSub,2}(1:l);
    plot(zeros(l,1),'--');
    hold all
    plot(diffPupil, 'color', 'k', 'linewidth', 1)
    title('high minus low');
    xlabel('time (ms)')
    xticks(0:1000:8000);
    xticklabels(2*[0:1000:8000]);
    ylabel('difference in pupil size (pixels)');
    set(gcf, 'position', [130 130 600 1000]);
end

%%

figure(6)
clf
rows=2;
cols=5;

subplot(rows,cols,1)
plot(allMeanPupil') %original mean pupil size
title('original mean pupil')
miny = min(allMeanPupil(:)) -400;
maxy = max(allMeanPupil(:)) + 400;
y = [miny maxy];
ylim(y);

subplot(rows,cols,cols+1)
plot(allResidMeanPupil') %mean pupil size after removing microsaccades convolved with pupil function
title('pupil after regressing out ms')
ylim(y);
% ylim([miny maxy]);


subplot(rows,cols, 2)
plot(allMeanPupil(1,:)'-allMeanPupil(2,:)');%difference between high and low
title('high minus low')
y2 = ylim;

subplot(rows,cols, cols+2)
plot(allResidMeanPupil(1,:)'-allResidMeanPupil(2,:)');%difference between high and low after removing ms
title('high minus low after regression')
ylim(y2);


subplot(rows,cols,rows*cols)
plot(allRegressMeanPupil') %what we regressed out
title('mean regressed out')




toc
