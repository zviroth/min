close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

nperms=1000;
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
end
zScoreString = '';
if toZscore
    zScoreString = '_zscored';
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
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
ROIs = [1:4];
% ROIs = [1 2];



%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;



%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd

% goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% goodSubs = [1 3 5:10 13]; %excluding 4 subjects with highest amplitude difference


clear smallSubAmp
for iRoi = 1:length(roiNames)
    for iSub=1:length(goodSubs)
        subTrialTC = [];
        for rwd=1:2
            subTrialTC = [subTrialTC subTrialResponse{goodSubs(iSub),iRoi,rwd}];
            
            numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);
        end
        subTrialAmp = std(subTrialTC);
        [sortedSubTrialAmp sortedInd] = sort(subTrialAmp,'descend');
        subTrialResponse{goodSubs(iSub),iRoi,1} = subTrialTC(:,sortedInd(1:numTrials(iSub,1)));
        subTrialResponse{goodSubs(iSub),iRoi,2} = subTrialTC(:,sortedInd(numTrials(iSub,1)+1:end));
        
        propTrueTrial(iRoi,iSub,1) = sum(sortedInd(1:numTrials(iSub,1))<=numTrials(iSub,1))/numTrials(iSub,1);
        propTrueTrial(iRoi,iSub,2) = sum(sortedInd(numTrials(iSub,1)+1:end)>numTrials(iSub,1))/numTrials(iSub,2);
        
    end
end
clear numTrials subTrialAmp


for iSub = 1:length(goodSubs)%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(roiNames)
        for rwd=1:2
%             reshapedTrials = reshape(subTrialResponse{goodSubs(iSub),iRoi,rwd},10,[]);
            reshapedTrials = subTrialResponse{goodSubs(iSub),iRoi,rwd};
           temp = fft(reshapedTrials);
           roiBinFftAmp = abs(temp(2,:));
           roiBinFftBaseline = abs(temp(1,:));
           
           %remove FFT component 1
%            reshapedTrials = reshapedTrials - repmat(mean(reshapedTrials),10,1);

           %remove FFT component 2
%            reshapedTrials = reshapedTrials./repmat(roiBinFftAmp,10,1);
%            reshapedTrials = (reshapedTrials - repmat(mean(reshapedTrials),10,1))./repmat(roiBinFftAmp,10,1) + repmat(mean(reshapedTrials),10,1);
           
           subTrialResponse{goodSubs(iSub),iRoi,rwd} = reshapedTrials;
    
        %trial-by-trial variability per subject
           subTimepointStd(iSub,iRoi,rwd,:) = std(reshapedTrials,0,2);
           temp = fft(reshapedTrials);
           roiBinFftAmp = abs(temp(2,:));
           roiBinFftPh = angle(temp(2,:));
           roiBinFftOther = sum(abs(temp([1 3:end],:)));
%            roiBinFftBaseline = abs(temp(1,:));
           roiBinFftBaseline = mean(reshapedTrials);
           subBinFftAmpStd(iSub,iRoi,rwd) = std(roiBinFftAmp);
           subBinFftPhStd(iSub,iRoi,rwd) = circ_std(roiBinFftPh');
           subBinFftOtherStd(iSub,iRoi,rwd) = std(roiBinFftOther);
           subBinFftBaselineStd(iSub,iRoi,rwd) = std(roiBinFftBaseline);

            if iSub==goodSubs(1)
                allTrials{iRoi,rwd} = reshapedTrials;
            else
                allTrials{iRoi,rwd} = [allTrials{iRoi,rwd} reshapedTrials];
            end
        end
     end
end


for rwd=1:2
    for iRoi = 1:length(roiNames)
        meanResponse(iRoi,rwd,:) = mean(subResponse(goodSubs,iRoi,rwd,:));
        stdResponse(iRoi,rwd,:) = std(subResponse(goodSubs,iRoi,rwd,:));
%         plot(squeeze(meanResponse(rwd,iRoi,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
%         hold on
        meanTimepointStd(iRoi, rwd,:) = mean(subTimepointStd(:,iRoi,rwd,:));%mean trial-by-trial variability
        stdTimepointStd(iRoi, rwd,:) = std(subTimepointStd(:,iRoi,rwd,:));%std of trial-by-trial variability
    end
end


%% PERMUTATIONS
tic

for iSub = 1:length(goodSubs)%length(subdirs)
    for rwd=1:2

        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!
        temp = trialResponse{goodSubs(iSub),rwd}(:,2:end-1);
%         temp = trialResponse{goodSubs(iSub),rwd}(:);
        trialResponseVec = temp(:);
%         temp = trialRT{goodSubs(iSub),rwd}(:);
        temp = trialRT{goodSubs(iSub),rwd}(:,2:end-1);
        trialRTvec{iSub,rwd} = temp(trialResponseVec>0);%only trials with a response, we're assuming onlyCorrect==0

%         if rwd==1
%             firstTrial(rwd)=1;
%         else %rwd==2
%             firstTrial(rwd)=numTrials(iSub,1)+1;
%         end
    end
    firstTrial(1)=1;
    firstTrial(2)=numTrials(iSub,1)+1;
    
    for iRoi=ROIs%length(roiNames)
        roiTrials{iRoi} = [subTrialResponse{goodSubs(iSub),iRoi,1} subTrialResponse{goodSubs(iSub),iRoi,2}];
    end
    subRT = [trialRTvec{iSub,1}; trialRTvec{iSub,2}];
    
    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));
        for iRoi= ROIs%length(roiNames)
            for rwd=1:2
                permTrials = roiTrials{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));% 10 timepoints X #trials
                permSubMeanTC(iSub,iRoi,rwd,p,:) = mean(permTrials,2);%average timecourse over trials
                permSubResp(iSub,iRoi,rwd,p) = std(squeeze(permSubMeanTC(iSub,iRoi,rwd,p,:)));%amplitude of task response
                
%                 temp = abs(fft(permSubMeanTC(iSub,iRoi,rwd,p,:)));
%                 permSubFFT(iSub,iRoi,rwd,p) = temp(ntrials+1);

                permSubStd(iSub,iRoi,rwd,p,:) = std(permTrials,0,2);%variability timecourse over trials
                
                %trial-by-trial FFT amp and phase variability per subject
                temp = fft(permTrials);
                permFftAmp = abs(temp(2,:));
                permFftPh = angle(temp(2,:));
                permFftOther = sum(abs(temp([1 3:end],:)));
                permSubAmpStd(iSub,iRoi,rwd,p) = std(permFftAmp);
                permSubPhStd(iSub,iRoi,rwd,p) = circ_std(permFftPh');
                permSubOtherStd(iSub,iRoi,rwd,p) = std(permFftOther);
                %                            permFftBaseline = abs(temp(1,:));
                permFftBaseline = mean(permTrials);
                permSubBaselineStd(iSub,iRoi,rwd,p) = std(permFftBaseline);
                

            end
        end
        for rwd=1:2
                %RT variability
                permRT = subRT(randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1),:);
                permSubRTvar(iSub,rwd,p) = std(permRT);
                permSubRT(iSub,rwd,p) = mean(permRT);
        end
    end
    
    for iRoi= ROIs%length(roiNames)
        for rwd = 1:2
            realSubMeanTC(iSub,iRoi,rwd,:) = mean(subTrialResponse{goodSubs(iSub),iRoi,rwd},2);
        end
        realSubDiff(iSub,iRoi) = std(realSubMeanTC(iSub,iRoi,1,:)) - std(realSubMeanTC(iSub,iRoi,2,:));%difference between high std and low std
        permSubDiff(iSub,iRoi,:) = squeeze(permSubResp(iSub,iRoi,1,:) - permSubResp(iSub,iRoi,2,:));%difference between high std and low std
        pVal_sub(iSub, iRoi) = sum(permSubDiff(iSub,iRoi,:) > realSubDiff(iSub,iRoi))/nperms;

        %FFT variability
        realSubAmpStdDiff(iSub,iRoi) = subBinFftAmpStd(iSub,iRoi,2) - subBinFftAmpStd(iSub,iRoi,1);
        realSubPhStdDiff(iSub,iRoi) = subBinFftPhStd(iSub,iRoi,2) - subBinFftPhStd(iSub,iRoi,1);
        realSubOtherStdDiff(iSub,iRoi) = subBinFftOtherStd(iSub,iRoi,2) - subBinFftOtherStd(iSub,iRoi,1);
        realSubBaselineStdDiff(iSub,iRoi) = subBinFftBaselineStd(iSub,iRoi,2) - subBinFftBaselineStd(iSub,iRoi,1);
        
        
        permSubAmpStdDiff(iSub,iRoi,:) = squeeze(permSubAmpStd(iSub,iRoi,2,:) - permSubAmpStd(iSub,iRoi,1,:));
        permSubPhStdDiff(iSub,iRoi,:) = squeeze(permSubPhStd(iSub,iRoi,2,:) - permSubPhStd(iSub,iRoi,1,:));
        permSubOtherStdDiff(iSub,iRoi,:) = squeeze(permSubOtherStd(iSub,iRoi,2,:) - permSubOtherStd(iSub,iRoi,1,:));
        permSubBaselineStdDiff(iSub,iRoi,:) = squeeze(permSubBaselineStd(iSub,iRoi,2,:) - permSubBaselineStd(iSub,iRoi,1,:));

    end
    for rwd=1:2
        %RT variability
        realSubRTvar(iSub,rwd) = std(trialRTvec{iSub,rwd});
        %RT
        realSubRT(iSub,rwd) = mean(trialRTvec{iSub,rwd});
    end
    
end
hypoth = pVal_sub<0.05;
% pVal_sub(:,ROIs);

meanRealAmpStdDiff = squeeze(mean(realSubAmpStdDiff));%mean over subjects
meanRealPhStdDiff = squeeze(mean(realSubPhStdDiff));%mean over subjects
meanRealOtherStdDiff = squeeze(mean(realSubOtherStdDiff));%mean over subjects
meanRealBaselineStdDiff = squeeze(mean(realSubBaselineStdDiff));%mean over subjects

meanPermAmpStdDiff = squeeze(mean(permSubAmpStdDiff));%mean over subjects
meanPermPhStdDiff = squeeze(mean(permSubPhStdDiff));%mean over subjects
meanPermOtherStdDiff = squeeze(mean(permSubOtherStdDiff));%mean over subjects
meanPermBaselineStdDiff = squeeze(mean(permSubBaselineStdDiff));%mean over subjects



for iRoi= ROIs
    %RFX
    meanPermDiff(iRoi,:) = mean(permSubDiff(:,iRoi,:));%average over subjects
    meanRealDiff(iRoi) = mean(realSubDiff(:,iRoi));%average over subjects
    pVal_rfx(iRoi) = sum(meanPermDiff(iRoi,:) > meanRealDiff(iRoi))/nperms;
    
%     %FFT RFX
%     meanPermFftDiff(iRoi,:) = mean(permSubFftDiff(:,iRoi,:));%average over subjects
%     meanRealFftDiff(iRoi) = mean(realSubFftDiff(:,iRoi));%average over subjects
%     pVal_fft_rfx(iRoi) = sum(meanPermFftDiff(iRoi,:) > meanRealFftDiff(iRoi))/nperms;
    
    %FFX
    for rwd = 1:2
        meanRealTC(iRoi,rwd,:) = squeeze(mean(realSubMeanTC(:,iRoi,rwd,:)));%average timecourse over subjects
        meanPermTC(iRoi,rwd,:,:) = squeeze(mean(permSubMeanTC(:,iRoi,rwd,:,:)));%average timecourse over subjects
    end
    meanRealDiffStd(iRoi) = std(squeeze(meanRealTC(iRoi,1,:))) - std(squeeze(meanRealTC(iRoi,2,:)));
    for p=1:nperms
        meanPermDiffStd(iRoi,p) = std(meanPermTC(iRoi,1,p,:)) - std(meanPermTC(iRoi,2,p,:));
    end
    pVal_ffx(iRoi) = sum(meanPermDiffStd(iRoi,:) > meanRealDiffStd(iRoi))/nperms;
    
    %variability
    permStd(iRoi,:,:,:) = mean(permSubStd(:,iRoi,:,:,:));%mean over subjects. permStd(iRoi,iRoi,rwd,p)
    for t=1:trialLength
       pVal_var_timecourse(iRoi,t) = sum( (permStd(iRoi,2,:,t)-permStd(iRoi,1,:,t)) > (meanTimepointStd(iRoi,2,t)-meanTimepointStd(iRoi,1,t)))/nperms;
    end
    pVal_var(iRoi) = sum( mean(permStd(iRoi,2,:,:)-permStd(iRoi,1,:,:),4) > mean(meanTimepointStd(iRoi,2,:)-meanTimepointStd(iRoi,1,:),3) )/nperms;
    
    %FFT variability
    pVal_ampVar(iRoi) = sum(meanPermAmpStdDiff(iRoi,:) > meanRealAmpStdDiff(iRoi))/nperms;
    pVal_phVar(iRoi) = sum(meanPermPhStdDiff(iRoi,:) > meanRealPhStdDiff(iRoi))/nperms;
    pVal_otherVar(iRoi) = sum(meanPermOtherStdDiff(iRoi,:) > meanRealOtherStdDiff(iRoi))/nperms;
    pVal_baselineVar(iRoi) = sum(meanPermBaselineStdDiff(iRoi,:) > meanRealBaselineStdDiff(iRoi))/nperms;

end
pVal_rfx
pVal_ffx
% pVal_fft_rfx
pVal_var
% pVal_var_timecourse*trialLength

pVal_ampVar
pVal_phVar
pVal_otherVar
pVal_baselineVar

%RT variability
permRTvar = squeeze(mean(permSubRTvar));
permRTvarDiff = permRTvar(2,:)-permRTvar(1,:);
realRTvar = mean(realSubRTvar);
realRTvarDiff = realRTvar(2) - realRTvar(1);
pVal_RTvar = sum(permRTvarDiff >= realRTvarDiff)/nperms;
%RT
permMeanRT = squeeze(mean(permSubRT));
permRTdiff = permMeanRT(2,:)-permMeanRT(1,:);
realMeanRT = mean(realSubRT);
realRTdiff = realMeanRT(2) - realMeanRT(1);
pVal_RT = sum(permRTdiff >= realRTdiff)/nperms;




%%
histbins=1000;
rows=length(roiNames);
cols = length(goodSubs);
% i=i+1;
% figure(i)
% clf
% for iSub = 1:length(goodSubs)
%     for iRoi= ROIs
%         subplot(rows,cols,iSub + (iRoi-1)*cols)
%         histogram(squeeze(permSubDiff(iSub,iRoi,:)),histbins); hold all
%         vline(prctile(permSubDiff(iSub,iRoi,:),95),'k');
%         vline(realSubDiff(iSub,iRoi),'r');
%     end
% end
        
        
%%
i=i+1;
figure(i)
clf
rows=length(roiNames);
cols = 3;
for iRoi= ROIs

    %phase
    subplot(rows,cols,1+(iRoi-1)*cols)
        histogram(squeeze(meanPermPhStdDiff(iRoi,:)),histbins); hold all
        vline(prctile(meanPermPhStdDiff(iRoi,:),95),'k');
        vline(meanRealPhStdDiff(iRoi),'r');
        title(['phase: ' roiNames{iRoi}]);
    %amplitude
        subplot(rows,cols,2+(iRoi-1)*cols)
        histogram(squeeze(meanPermAmpStdDiff(iRoi,:)),histbins); hold all
        vline(prctile(meanPermAmpStdDiff(iRoi,:),95),'k');
        vline(meanRealAmpStdDiff(iRoi),'r');
        title(['amplitude: ' roiNames{iRoi}]);
    %other FFT
        subplot(rows,cols,3+(iRoi-1)*cols)
        histogram(squeeze(meanPermOtherStdDiff(iRoi,:)),histbins); hold all
        vline(prctile(meanPermOtherStdDiff(iRoi,:),95),'k');
        vline(meanRealOtherStdDiff(iRoi),'r');
        title(['other FFT: ' roiNames{iRoi}]);
end



%%
iRoi=2;
i=i+1;
figure(i)
clf
cols = length(goodSubs);
rows=2;
for iSub=1:length(goodSubs)
    for rwd=1:2
        subplot(rows,cols,(rwd-1)*cols+iSub)
        plot(subTrialResponse{goodSubs(iSub),iRoi,rwd});
        hold on
        plot(mean(subTrialResponse{goodSubs(iSub),iRoi,rwd},2),'k','linewidth',2);
    end
end
%%
i=i+1;
figure(i)
clf
cols = length(goodSubs);
rows=1;
for iSub=1:length(goodSubs)
    subplot(rows,cols,iSub)
    for rwd=1:2
        
%         plot(subTrialResponse{goodSubs(iSub),iRoi,rwd});
%         hold on
        plot(mean(subTrialResponse{goodSubs(iSub),iRoi,rwd},2),'color',plotColors{rwd},'linewidth',2);
%         plot(squeeze(subResponse(goodSubs(iSub),iRoi,rwd,:)), 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',20);

        hold on
    end
end


%%
globalMean{iSub,rwd}\roiTC{iSub,iRoi,rwd}.tSeries';

% %%
% i=i+1;
% figure(i)
% clf
% rows=1;
% cols=length(roiNames);
% for iRoi = 1:length(roiNames)
%     subplot(rows,cols,iRoi)
%     for rwd=1:2
%         dsErrorsurface(1:trialLength, squeeze(meanStd(iRoi,rwd,:)), squeeze(stdStd(iRoi,rwd,:))./sqrt(size(subStd,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
%         hold on
%         plot(1:10, squeeze(meanStd(iRoi,rwd,:)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
%     end
% end


clear smallSubAmp
for iRoi = 1:length(roiNames)
    for iSub=1:length(goodSubs)
        for rwd=1:2
            smallSubAmp(iSub,iRoi,rwd) = std(squeeze(subResponse(goodSubs(iSub),iRoi,rwd,:)));
        end
    end
end


iRoi=2;

smallSubDiff = smallSubAmp(:,:,1) - smallSubAmp(:,:,2);
subStdDiff = mean(subTimepointStd(:,:,2,:),4) - mean(subTimepointStd(:,:,1,:),4);
[r p] = corr(realSubAmpStdDiff(:,iRoi),subStdDiff(:,iRoi))
[r p] = corr(realSubPhStdDiff(:,iRoi),subStdDiff(:,iRoi))
[r p] = corr(realSubPhStdDiff(:,iRoi),realSubAmpStdDiff(:,iRoi))

[r p] = corr(smallSubDiff(:,iRoi),subStdDiff(:,iRoi))
[r p] = corr(smallSubDiff(:,iRoi),realSubAmpStdDiff(:,iRoi))
[r p] = corr(smallSubDiff(:,iRoi),realSubPhStdDiff(:,iRoi))




%%  correlate amplitude difference and phase variability with RT variability
% saveFolder= '/Volumes/MH02086153MACDT-Drobo/allMinBehavioral/';
% load([saveFolder 'RTvariability' onlyCorrectString '.mat'], 'realRTvarDiff','permRTvarDiff','pVal_RTvar',...
%     'realRTvar','realSubRTvar','trialRT','permRTvar','permSubRTvar','subFolders','nperms',...
%     'pVal_std_rfx','pVal_baseline_rfx','pVal_var','meanRealStdDiff','realSubStdDiff','meanPermVarDiff',...
%     'permBaselineDiff','permSubVar','permSubStdDiff','realBaselineDiff','realStdDiff');
pVal_RTvar
[r p] = corr(realSubRTvar(:,2)-realSubRTvar(:,1),realSubPhStdDiff(:,iRoi))
[r p] = corr(realSubRTvar(:,2)-realSubRTvar(:,1),smallSubDiff(:,iRoi))

pVal_RT
[r p] = corr(realSubRT(:,2)-realSubRT(:,1),realSubPhStdDiff(:,iRoi))
[r p] = corr(realSubRT(:,2)-realSubRT(:,1),smallSubDiff(:,iRoi))






toc

