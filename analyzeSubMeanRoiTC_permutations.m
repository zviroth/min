dataFolder = '/Volumes/MH02086153MACDT-Drobo/decodingAnalysis/rwd/';
onlyCorrect=1;
onlyCorrectString = '';
if onlyCorrect
    onlyCorrectString = '_correct';
end
load([dataFolder 'subMeanRoiTC' onlyCorrectString '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', 'meanResponse', 'stdResponse',...
    'roiTC', 'allTrials', ...
    'subdirs', 'roiNames','subTrialResponse','subRunResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'rwdPupil','meanPupil','expName');

clear subMeanResponse trialStd meanTrialStd runStd meanRunStd subRwdStd trialFFTamp trialFFTphase meanTrialFFTamp meanTrialFFTphase 
clear stdTrialFFTphase runFFTamp runFFTphase meanRunFFTamp meanRunFFTphase stdRunFFTphase subFFTamp subFFTphase
clear trialMinMax meanTrialMinMax runMinMax meanRunMinMax subRwdMinMax legendCellArray
clear stdRunMax stdTrialMax
clear groupLabels pVal
clear yError yMin yMax


% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
ROIs = [1 2 7 8];
ROIs = [1 2];

% ROIs = 1:4;

for i=1:length(ROIs)
    groupLabels{i} = roiNames{ROIs(i)};
end
legendHighLow = [];
for iRoi= 1:length(groupLabels)
    legendHighLow{(iRoi-1)*2+1} = [groupLabels{iRoi} ': H'];
    legendHighLow{iRoi*2} = [groupLabels{iRoi} ': L'];
end
% ROIs = [1:2];

% roiNames(5:end) = [];
% roiNames([1:2 4:end]) = [];
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

for iSub = 1:length(subdirs)
    subplot(rows,cols,iSub)
    for iRoi= ROIs%1:length(roiNames)
        
        for rwd=1:2
            
%             STD
            trialStd{iSub,iRoi,rwd} = std(subTrialResponse{iSub,iRoi,rwd});%std per trial
            meanTrialStd(iSub,iRoi,rwd) = mean(trialStd{iSub,iRoi,rwd});
            runStd{iSub,iRoi,rwd} = std(subRunResponse{iSub,iRoi,rwd});%std per run
            meanRunStd(iSub,iRoi,rwd) = mean(runStd{iSub,iRoi,rwd});
            subRwdStd(iSub,iRoi,rwd) = std(squeeze(subResponse(iSub,iRoi,rwd,:)));
            
            [maxValTrial{iSub,iRoi,rwd} maxIndTrial{iSub,iRoi,rwd}] = max(subTrialResponse{iSub,iRoi,rwd});
            [maxValRun{iSub,iRoi,rwd} maxIndRun{iSub,iRoi,rwd}] = max(subRunResponse{iSub,iRoi,rwd});
            
            [s s0] =  circ_std(maxIndTrial{iSub,iRoi,rwd}'*2*pi/10);%radians
            stdTrialMax(iSub,iRoi,rwd) = s0;
            [s s0] =  circ_std(maxIndRun{iSub,iRoi,rwd}'*2*pi/10);%radians
            stdRunMax(iSub,iRoi,rwd) = s0;
            
%             FFT
            trialFFT = fft(subTrialResponse{iSub,iRoi,rwd});
            trialFFTamp{iSub,iRoi,rwd} = abs(trialFFT(2,:));
            trialFFTphase{iSub,iRoi,rwd} = angle(trialFFT(2,:));
            meanTrialFFTamp(iSub,iRoi,rwd) = mean(trialFFTamp{iSub,iRoi,rwd});
%             meanTrialFFTphase(iSub,iRoi,rwd) = mean(trialFFTphase{iSub,iRoi,rwd});
%             stdTrialFFTphase(iSub,iRoi,rwd) = std(trialFFTphase{iSub,iRoi,rwd});
            [s s0] =  circ_std(trialFFTphase{iSub,iRoi,rwd}');
            stdTrialFFTphase(iSub,iRoi,rwd) = s0;
            stdTrialFFTamp(iSub,iRoi,rwd) = std(trialFFTamp{iSub,iRoi,rwd},0,2);
            
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
            title([subdirs(iSub).name(1:5) newline expName{iSub,rwd}]);
            hold on
        end
        %FFT of mean response, across reward type
        subMeanResponse{iSub,iRoi} = mean(squeeze(subResponse(iSub,iRoi,:,:)));%averaged over reward type
        temp = fft(subMeanResponse{iSub,iRoi});
        subMeanFFT(iSub,iRoi) = temp(2);
        subMeanFFTamp(iSub,iRoi) = abs(temp(2));
        subMeanFFTph(iSub,iRoi) = angle(temp(2));
        
    end
%     %correlation between V1-DMN amplitudes 
%     for rwd=1:2
%         trialStdCorr{iSub} = corr(trialStd{iSub,1,rwd}', trialStd{iSub,2,rwd}');
%         runStdCorr{iSub} = corr(runStd{iSub,1,rwd}', runStd{iSub,2,rwd}');
%     end
end
legend(legendHighLow);
set(gcf,'position',[250 500 1200 600]);

%% DIFFERENT MEASURES OF RESPONSE, PER SUBJECT/ROI
figure(2)
clf
rows=3;
cols=3;
iRoi=1;
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
plot(phaseNorm*(stdTrialFFTphase(:,iRoi,1)-stdTrialFFTphase(:,iRoi,2)),'b.', 'linewidth', 1);
% corr(meanTrialFFTamp(:,1,1)-meanTrialFFTamp(:,1,2), stdTrialFFTphase(:,1,2)-stdTrialFFTphase(:,1,1))
subplot(rows,cols,5)
plot(meanRunFFTamp(:,:,1)-meanRunFFTamp(:,:,2), 'linewidth', 1); title('runs FFT');
hold on
plot(phaseNorm*(stdRunFFTphase(:,iRoi,1)-stdRunFFTphase(:,iRoi,2)),'b.', 'linewidth', 1);
% corr(meanRunFFTamp(:,1,1)-meanRunFFTamp(:,1,2),stdRunFFTphase(:,1,2)-stdRunFFTphase(:,1,1))
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
legend(groupLabels);
set(gcf,'position',[200 150 900 700]);
%% BAR PLOTS
clear meanRwdStd
figure(3)
clf
xshift=0.2;
scatterSize = 30;
linewidth = 3;
linelength = 0.2;
%line per subject per ROI, connecting low and high reward
for i= 1:length(ROIs)%1:length(roiNames)
    iRoi = ROIs(i);
    line([i-0.5*xshift i+0.5*xshift], [subRwdStd(:,iRoi,1) subRwdStd(:,iRoi,2)],'linewidth',linewidth/2);
    hold all
    pval(iRoi) = ranksum(subRwdStd(:,iRoi,1), subRwdStd(:,iRoi,2));
end
%black bars for average across subjects
for rwd=1:2
    meanRwdStd(:,rwd) = mean(subRwdStd(:,:,rwd)); %per ROI
%     stdRwdStd(:,rwd) = std(subRwdStd(:,:,rwd)); %per ROI
    for i= 1:length(ROIs)%iRoi= ROIs%1:length(roiNames)
        iRoi = ROIs(i);
        scatterCenter = i+(rwd-1.5)*xshift;
        scatter(scatterCenter*ones(length(subdirs),1), subRwdStd(:,iRoi,rwd),scatterSize, plotColors{rwd});
        hold all
        line([scatterCenter-linelength/2 scatterCenter+linelength/2], [meanRwdStd(iRoi,rwd) meanRwdStd(iRoi,rwd)], 'color','k','linewidth',linewidth);
    end 
end
%black line connecting between means
for i= 1:length(ROIs)%iRoi= ROIs%1:length(roiNames)
    iRoi = ROIs(i);
    line([i-0.5*xshift i+0.5*xshift], [meanRwdStd(iRoi,1) meanRwdStd(iRoi,2)],'linewidth',linewidth/2, 'color', 'k','linewidth',linewidth);
end


xlim([0 length(groupLabels)+1]);
xticks(1:length(groupLabels));
xticklabels(groupLabels);



%% average timecourse for all ROIs
figure(5)
clf
cols=ceil(length(subdirs)/2);
rows = ceil(length(subdirs)/cols);
for iRoi= ROIs%1:length(roiNames)
    for rwd=1:2
        plot(squeeze(meanResponse(rwd,iRoi,:)),plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
        hold on
    end
end
legend(legendHighLow);

%% plot all run average timecourses for right V1
figure(6)
rwdColors = {'b','r'};
roiNames = {'rV1_eccen8', 'lV1_eccen8','rV2_eccen8','lV2_eccen8','rV3_eccen8','lV3_eccen8','rh_17Networks_16','lh_17Networks_16'};
rV1roi = 0;
for iRoi=1:length(roiNames)
    if strcmp('rV1_eccen8', roiNames{iRoi})
        rV1roi = iRoi;
    end
end
for iSub = 1:length(subdirs)
    subplot(rows,cols,iSub)
    for rwd=1:2
        plot(subRunResponse{iSub,rV1roi,rwd},'linestyle','-','marker','.','color',rwdColors{rwd});
        hold on
    end
end



%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials
% subject 4 = 22
% subject 16 = outlier of RT and RT std
goodSubs = [1:3 5:18];
goodSubs = [1:7 9:10 12 15];
goodSubs = [8 11 13:14 16:17];
goodSubs = 1:18;
for iSub = goodSubs%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(roiNames)
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

nperms=10000;
for rwd=1:2
    numTrials(rwd) = size(allTrials{1,rwd},2);
    if rwd==1
        firstTrial(rwd)=1;
    else %rwd==2
        firstTrial(rwd)=numTrials(1)+1;
    end
end
for iRoi=ROIs%length(roiNames)
   roiTrials{iRoi} = [allTrials{iRoi,1} allTrials{iRoi,2}]; 
end

for p=1:nperms
    randOrder = randperm(numTrials(1)+numTrials(2));
    for iRoi= ROIs%length(roiNames)
        for rwd=1:2
            permTrials = roiTrials{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(rwd)-1));
            permResp(iRoi,rwd,p) = std(mean(permTrials,2));
        end
    end
end
for iRoi= ROIs%length(roiNames)
    realDiff(iRoi) = std(mean(allTrials{iRoi,1},2)) - std(mean(allTrials{iRoi,2},2)); 
%     realDiff(iRoi) = meanRwdStd(iRoi,1) - meanRwdStd(iRoi,2);
    permDiff = squeeze(permResp(iRoi,1,:) - permResp(iRoi,2,:));
    pVal(iRoi) = sum(permDiff>realDiff(iRoi))/nperms;
end
    
pVal(ROIs)

for iRoi = ROIs
   [h pVal_ttest(iRoi)] = ttest(subRwdStd(:,iRoi,1),subRwdStd(:,iRoi,2));
end


%% SINGLE SUBJECT PERMUTATIONS
tic
pVal_sub = ones(length(subdirs),2);
for iSub = goodSubs%length(subdirs)
    clear subTrials permStd
    for rwd=1:2
%         for iRoi=1:length(roiNames)
%             subTrials{iRoi,rwd} = reshape(subTrialResponse{iSub,iRoi,rwd},10,[]);
%         end
        numTrials(iSub,rwd) = size(subTrialResponse{iSub,1,rwd},2);
        if rwd==1
            firstTrial(rwd)=1;
        else %rwd==2
            firstTrial(rwd)=numTrials(iSub,1)+1;
        end
    end
    for iRoi=ROIs%length(roiNames)
        roiTrials{iRoi} = [subTrialResponse{iSub,iRoi,1} subTrialResponse{iSub,iRoi,2}];
        % FFT of each trial
        roiTrialsFFT{iRoi} = fft(roiTrials{iRoi});
    end
    
    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));
        for iRoi= ROIs%length(roiNames)
            for rwd=1:2
                permTrials = roiTrials{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));
                permSubResp(iSub,iRoi,rwd,p,:) = mean(permTrials,2);
                permResp(iSub,iRoi,rwd,p) = std(permSubResp(iSub,iRoi,rwd,p,:),0,5);
                permFFT = roiTrialsFFT{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));
%                 permFFT = fft(permTrials);
                permFFTamp = abs(permFFT(2,:));
                permFFTphase = angle(permFFT(2,:));
                [s s0] =  circ_std(permFFTphase');
                stdPermFFTphase(iRoi,rwd,p) = s0;
                stdPermFFTamp(iRoi,rwd,p) = std(permFFTamp);

            end
        end
    end
    
    for iRoi= ROIs%length(roiNames)
        for rwd = 1:2
            realSubResp(iSub,iRoi,rwd,:) = mean(subTrialResponse{iSub,iRoi,rwd},2);
        end
        realDiff(iRoi) = std(realSubResp(iSub,iRoi,2,:)) - std(realSubResp(iSub,iRoi,2,:));
        %     realDiff(iRoi) = meanRwdStd(iRoi,1) - meanRwdStd(iRoi,2);
        permDiff = squeeze(permResp(iSub,iRoi,1,:) - permResp(iSub,iRoi,2,:));
        pVal_sub(iSub,iRoi) = sum(permDiff>realDiff(iRoi))/nperms;
        
        realFFTphaseStdDiff(iSub,iRoi) = stdTrialFFTphase(iSub,iRoi,1) - stdTrialFFTphase(iSub,iRoi,2);
        permFFTphaseStdDiff(iSub,iRoi,:) = stdPermFFTphase(iRoi,1,:) - stdPermFFTphase(iRoi,2,:);
%         pVal_FFTphase(iSub,iRoi) = sum(permFFTphaseStdDiff < realFFTphaseStdDiff(iRoi))/nperms;
        
        realFFTampStdDiff(iSub,iRoi) = stdTrialFFTamp(iSub,iRoi,1) - stdTrialFFTamp(iSub,iRoi,2);
        permFFTampStdDiff(iSub,iRoi,:) = stdPermFFTamp(iRoi,1,:) - stdPermFFTamp(iRoi,2,:);
%         pVal_FFTamp(iSub,iRoi) = sum(permFFTampStdDiff < realFFTampStdDiff(iRoi))/nperms;
    end
end
pVal_sub(:,ROIs);
%%
%now combine across subjects for real and permuted data
for iRoi= ROIs%length(roiNames)
    meanPermFFTphaseStdDiff = mean(permFFTphaseStdDiff(:,iRoi,:));%mean across subjects
    meanRealFFTphaseStdDiff = mean(realFFTphaseStdDiff(:,iRoi));
    meanPermFFTampStdDiff = mean(permFFTampStdDiff(:,iRoi,:));%mean across subjects
    meanRealFFTampStdDiff = mean(realFFTampStdDiff(:,iRoi));
    
    
    pVal_FFTphase(iRoi) = sum(meanPermFFTphaseStdDiff < meanRealFFTphaseStdDiff)/nperms;
    pVal_FFTamp(iRoi) = sum(meanPermFFTampStdDiff < meanRealFFTampStdDiff)/nperms;
    
    for rwd=1:2
        for p=1:nperms
            meanPermResp(iRoi,rwd,p,:) = mean(squeeze(permSubResp(:,iRoi,rwd,p,:)) .* repmat(numTrials(:,rwd),1,10)); %averaging across subjects, weighted by #trials
            meanPermAmp(iRoi,rwd,p) = std(meanPermResp(iRoi,rwd,p,:));
            
        end
        meanRealResp(iRoi,rwd,:) = mean(squeeze(realSubResp(:,iRoi,rwd,:)) .* repmat(numTrials(:,rwd),1,10)); %averaging across subjects, weighted by #trials
        meanRealAmp(iRoi,rwd) = std(meanRealResp(iRoi,rwd,:));
        
    end
    meanRealDiff(iRoi) = meanRealAmp(iRoi,1) - meanRealAmp(iRoi,2);
    meanPermDiff(iRoi,:) = meanPermAmp(iRoi,1,:) - meanPermAmp(iRoi,2,:);
    pVal_subPerm(iRoi) = sum(meanPermDiff(iRoi,:) > meanRealDiff(iRoi))/nperms;
end
pVal_FFTphase
pVal_FFTamp
pVal_subPerm
toc
%%
figure(7)
clf
% meanRwdStd(roi,rwd)

% groupLabels = roiNames;
for rwd=1:2
    yError(:,rwd) = std(subRwdStd(:,ROIs,rwd));
    yMin(:,rwd) = min(subRwdStd(:,ROIs,rwd));
    yMax(:,rwd) = max(subRwdStd(:,ROIs,rwd));
end
yAxisMin = min(subRwdStd(:));
yAxisMax = max(subRwdStd(:));
dispValues=0;
mybar(meanRwdStd(ROIs,:), 'groupLabels',groupLabels, 'withinGroupLabels',{'H','L'} ,'yAxisMin', yAxisMin, 'yAxisMax', yAxisMax,'dispValues',dispValues,...
    'yError', yError, 'yMin', yMin, 'yMax', yMax);

%%
subMeanPropCorrect = cellfun(@mean, propCorrect);
[a pVal_propCorrect] = ttest(subMeanPropCorrect(:,1), subMeanPropCorrect(:,2));

%%
figure(8)
clf
minVal = min(min(real(subMeanFFT(:))),min(imag(subMeanFFT(:))));
maxVal = max(max(real(subMeanFFT(:))),max(imag(subMeanFFT(:))));
maxVal = max(abs(subMeanFFT(:)));
for i= 1:length(ROIs)%1:length(roiNames)
    subplot(1,length(ROIs),i);
    iRoi = ROIs(i);
    %     plot(real(subMeanFFT(:,iRoi)), imag(subMeanFFT(:,iRoi)),'.');
    
    %     hold all
%     axis([minVal maxVal minVal maxVal]);
    axis([-maxVal maxVal -maxVal maxVal]);
    line([-maxVal maxVal], [0  0],'color',[0 0 0],'linewidth', 2);
    line([0  0], [-maxVal maxVal],'color',[0 0 0],'linewidth', 2);
    hold all
    scatter(real(subMeanFFT(:,iRoi)), imag(subMeanFFT(:,iRoi)),50,'filled');
    axis image
end
%%
figure(9)
clf
for iSub = 1:length(subdirs)-1
    subplot(rows,cols,iSub)
    for rwd=1:2
        meanPupil{iSub,rwd} = nanmean(rwdPupil{iSub,rwd})';
        plot(meanPupil{iSub,rwd}, 'Color', plotColors{rwd}, 'linewidth', 1)
        hold on
    end
end
legend('high', 'low')

%% Mean proportion correct
figure(10)
clf
for iSub = 1:length(subdirs)
   for rwd=1:2
       meanPropCorrect(iSub,rwd) = mean(propCorrect{iSub,rwd});
       meanRT(iSub,rwd) = mean(trialRT{iSub,rwd,r});
       stdRT(iSub,rwd) = std(trialRT{iSub,rwd,r});
   end 
end
subPropDiff = meanPropCorrect(:,1) - meanPropCorrect(:,2);
subRTdiff = meanRT(:,1) - meanRT(:,2);
subRTstdDiff = stdRT(:,1) - stdRT(:,2);
plot(zscore(subPropDiff))
hold all
plot(zscore(subRTdiff))
plot(zscore(subRTstdDiff))
legend('prop correct', 'RT mean','RT std')

plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
% subRwdDiff = subRwdStd(:,:,1)-subRwdStd(:,:,2);

