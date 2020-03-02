close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

% randSeed = rng;
load(fullfile(dataFolder,'randSeed.mat'),'randSeed');
rng(randSeed);
% save(fullfile(dataFolder,'randSeed.mat'),'randSeed');


nperms=10000;
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
elseif onlyCorrect==0
    onlyCorrectString = '_validresponse';
end
zScoreString = '';
if toZscore
    zScoreString = '_zscored';
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
    'globalMean','regressBetasGlobal','runRwd',...
    'subRoiRuns','runMeanFFT',...
    'allVoxTrialResponse','allVoxTaskPhase','allVoxTaskAmp','allVoxTaskCo',...
    'voxTrials','voxGoodTrials','meanVoxTrial',...
    'maxRT');

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

goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% goodSubs = [1 3 5:10 13]; %excluding 4 subjects with highest amplitude difference

% ROIs = [1:4];
% ROIs = [1 2];
rwdString = {'H','L'};
%% ANALYZE GLOBAL MEAN SIGNAL
for iSub = 1:length(goodSubs)
    for rwd=1:2
         globalTrials = reshape(globalMean{iSub,rwd},trialLength,[]);
         meanGlobalTrial(iSub,rwd,:) = mean(globalTrials');
         varGlobalTrial(iSub,rwd,:) = std(globalTrials');
    end
end

for rwd=1:2
    plot(squeeze(mean(meanGlobalTrial(:,rwd,:))),'color',plotColors{rwd});
    hold on
end





i=0;
%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;

%% mean + STD of thresholds
subThresh = mean(subMeanThresh,2);%mean over rwd, which is already averaged over staircases and runs
mean(subThresh);
std(subThresh);


%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd

for iSub = 1:length(goodSubs)%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(roiNames)
        for rwd=1:2
%             reshapedTrials = reshape(subTrialResponse{goodSubs(iSub),iRoi,rwd},10,[]);
            reshapedTrials = subTrialResponse{goodSubs(iSub),iRoi,rwd};
           temp = fft(reshapedTrials);
           roiFftAmp = abs(temp(2,:));
           roiFftBaseline = abs(temp(1,:));
           
           %remove FFT component 1
%            reshapedTrials = reshapedTrials - repmat(mean(reshapedTrials),10,1);

           %remove FFT component 2
%            reshapedTrials = reshapedTrials./repmat(roiBinFftAmp,10,1);
%            reshapedTrials = (reshapedTrials - repmat(mean(reshapedTrials),10,1))./repmat(roiFftAmp,10,1) + repmat(mean(reshapedTrials),10,1);
           
           subTrialResponse{goodSubs(iSub),iRoi,rwd} = reshapedTrials;
    
        %trial-by-trial variability per subject
           subTimepointStd(iSub,iRoi,rwd,:) = std(reshapedTrials,0,2);
           temp = fft(reshapedTrials);
           roiFftAmp = abs(temp(2,:));
           roiFftPh = angle(temp(2,:));
           trialRoiFftPhVec{iSub,iRoi,rwd} = roiFftPh;
           trialRoiFftAmpVec{iSub,iRoi,rwd} = roiFftAmp;
           trialRoiAmpVec{iSub,iRoi,rwd} = std(reshapedTrials);
           roiFftOther = sum(abs(temp([1 3:end],:)));
%            roiBinFftBaseline = abs(temp(1,:));
           roiFftBaseline = mean(reshapedTrials);
           subFftAmpStd(iSub,iRoi,rwd) = std(roiFftAmp);
           subFftMeanAmp(iSub,iRoi,rwd) = mean(roiFftAmp);
           subFftMeanPh(iSub,iRoi,rwd) = circ_mean(roiFftPh');
           subMeanStd(iSub,iRoi,rwd) = mean(std(reshapedTrials));
           subStdVar(iSub,iRoi,rwd) = std(std(reshapedTrials));%variability of amplitude measured by std
           subMedianStd(iSub,iRoi,rwd) = median(std(reshapedTrials));
           
           subFftPhStd(iSub,iRoi,rwd) = circ_std(roiFftPh');
           subFftOtherStd(iSub,iRoi,rwd) = std(roiFftOther);
           subFftBaselineStd(iSub,iRoi,rwd) = std(roiFftBaseline);

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
        temp = trialCorrectness{goodSubs(iSub),rwd}(:,2:end-1);%exclude the first trial and the last trial
        trialCorrectVec{iSub,rwd} = temp(:); %include trials with no response
        subPropCorrect(iSub,rwd) = sum(trialCorrectVec{iSub,rwd})/length(trialCorrectVec{iSub,rwd});
        
        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!
        temp = trialResponse{goodSubs(iSub),rwd}(:,2:end-1);
%         temp = trialResponse{goodSubs(iSub),rwd}(:);
        trialResponseVec = temp(:);
%         temp = trialRT{goodSubs(iSub),rwd}(:);
        temp = trialRT{goodSubs(iSub),rwd}(:,2:end-1);
        allTrialRTvec = temp(:);
        
        if onlyCorrect==0
            trialRTvec{iSub,rwd} = temp(trialResponseVec>0 & allTrialRTvec>0 & allTrialRTvec<maxRT);%only trials with a response, we're assuming onlyCorrect==0
        elseif onlyCorrect==1 %only correct
            trialRTvec{iSub,rwd} = temp(trialCorrectVec{iSub,rwd}==1 & allTrialRTvec>0 & allTrialRTvec<maxRT);
        elseif onlyCorrect==2 %only incorrect
            trialRTvec{iSub,rwd} = temp(trialCorrectVec{iSub,rwd}==0 & trialResponseVec>0 & allTrialRTvec>0 & allTrialRTvec<maxRT);
        else
            trialRTvec{iSub,rwd} = temp(:);
        end

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
                
                %latency
                temp=fft(permSubMeanTC(iSub,iRoi,rwd,p,:));
                permSubMeanFftAmp(iSub,iRoi,rwd,p) = abs(temp(2));
                permSubMeanFftPh(iSub,iRoi,rwd,p) = angle(temp(2));
                
                
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
                
                permSubMeanAmp(iSub,iRoi,rwd,p) = mean(permFftAmp);
                permSubMeanPh(iSub,iRoi,rwd,p) = circ_mean(permFftPh');
                permSubMeanStd(iSub,iRoi,rwd,p) = mean(std(permTrials));
                permSubMedianStd(iSub,iRoi,rwd,p) = median(std(permTrials));
                
                permSubStdVar(iSub,iRoi,rwd,p) = std(std(permTrials));
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
        %difference in amplitude
        for rwd = 1:2
            realSubMeanTC(iSub,iRoi,rwd,:) = mean(subTrialResponse{goodSubs(iSub),iRoi,rwd},2);
            temp = fft(squeeze(realSubMeanTC(iSub,iRoi,rwd,:)));
            realSubMeanFftAmp(iSub,iRoi,rwd) = abs(temp(2));
            realSubMeanFftPh(iSub,iRoi,rwd) = angle(temp(2));
        end
        realSubDiff(iSub,iRoi) = std(realSubMeanTC(iSub,iRoi,1,:)) - std(realSubMeanTC(iSub,iRoi,2,:));%difference between high std and low std
        permSubDiff(iSub,iRoi,:) = squeeze(permSubResp(iSub,iRoi,1,:) - permSubResp(iSub,iRoi,2,:));%difference between high std and low std
        pVal_sub(iSub, iRoi) = sum(permSubDiff(iSub,iRoi,:) >= realSubDiff(iSub,iRoi))/nperms;
        
        %difference in latency
        realSubFftAmpDiff(iSub,iRoi) = realSubMeanFftAmp(iSub,iRoi,1) - realSubMeanFftAmp(iSub,iRoi,2);
        realSubFftPhDiff(iSub,iRoi) = circ_dist(realSubMeanFftPh(iSub,iRoi,2), realSubMeanFftPh(iSub,iRoi,1));
        permSubFftAmpDiff(iSub,iRoi,:) = permSubMeanFftAmp(iSub,iRoi,1,:) - permSubMeanFftAmp(iSub,iRoi,2,:);
        permSubFftPhDiff(iSub,iRoi,:) = circ_dist(permSubMeanFftPh(iSub,iRoi,2,:),permSubMeanFftPh(iSub,iRoi,1,:));
        pVal_fftPh_sub(iSub, iRoi) = sum(permSubFftPhDiff(iSub,iRoi,:) >= realSubFftPhDiff(iSub,iRoi))/nperms;
        
        %FFT variability
        realSubAmpStdDiff(iSub,iRoi) = subFftAmpStd(iSub,iRoi,2) - subFftAmpStd(iSub,iRoi,1);
        realSubPhStdDiff(iSub,iRoi) = subFftPhStd(iSub,iRoi,2) - subFftPhStd(iSub,iRoi,1);
        realSubOtherStdDiff(iSub,iRoi) = subFftOtherStd(iSub,iRoi,2) - subFftOtherStd(iSub,iRoi,1);
        realSubBaselineStdDiff(iSub,iRoi) = subFftBaselineStd(iSub,iRoi,2) - subFftBaselineStd(iSub,iRoi,1);
        
        realSubMeanFftAmpDiff(iSub,iRoi) = subFftMeanAmp(iSub,iRoi,2) - subFftMeanAmp(iSub,iRoi,1);
        realSubMeanFftPhDiff(iSub,iRoi) = squeeze(circ_dist(subFftMeanPh(iSub,iRoi,2),subFftMeanPh(iSub,iRoi,1)));
        realSubMeanStdDiff(iSub,iRoi) = subMeanStd(iSub,iRoi,2) - subMeanStd(iSub,iRoi,1);
        realSubMedianStdDiff(iSub,iRoi) = subMedianStd(iSub,iRoi,2) - subMedianStd(iSub,iRoi,1);
        
        realSubStdVarDiff(iSub,iRoi) = subStdVar(iSub,iRoi,2) - subStdVar(iSub,iRoi,1);
        
        permSubAmpStdDiff(iSub,iRoi,:) = squeeze(permSubAmpStd(iSub,iRoi,2,:) - permSubAmpStd(iSub,iRoi,1,:));
        permSubPhStdDiff(iSub,iRoi,:) = squeeze(permSubPhStd(iSub,iRoi,2,:) - permSubPhStd(iSub,iRoi,1,:));
        permSubOtherStdDiff(iSub,iRoi,:) = squeeze(permSubOtherStd(iSub,iRoi,2,:) - permSubOtherStd(iSub,iRoi,1,:));
        permSubBaselineStdDiff(iSub,iRoi,:) = squeeze(permSubBaselineStd(iSub,iRoi,2,:) - permSubBaselineStd(iSub,iRoi,1,:));
        
        permSubMeanAmpDiff(iSub,iRoi,:) = squeeze(permSubMeanAmp(iSub,iRoi,2,:) - permSubMeanAmp(iSub,iRoi,1,:));
        permSubMeanPhDiff(iSub,iRoi,:) = squeeze(circ_dist(permSubMeanPh(iSub,iRoi,2,:),permSubMeanPh(iSub,iRoi,1,:)));
        permSubMeanStdDiff(iSub,iRoi,:) = squeeze(permSubMeanStd(iSub,iRoi,2,:) - permSubMeanStd(iSub,iRoi,1,:));
        permSubMedianStdDiff(iSub,iRoi,:) = squeeze(permSubMedianStd(iSub,iRoi,2,:) - permSubMedianStd(iSub,iRoi,1,:));
        
        permSubStdVarDiff(iSub,iRoi,:) = squeeze(permSubStdVar(iSub,iRoi,2,:) - permSubStdVar(iSub,iRoi,1,:));
    end
    for rwd=1:2
        %RT variability
        realSubRTvar(iSub,rwd) = std(trialRTvec{iSub,rwd});
        %RT
        realSubRT(iSub,rwd) = mean(trialRTvec{iSub,rwd});
    end
    
end





%single subject p-values for timepoint variability
realSubTimepointVar = squeeze(mean(subTimepointStd,4));
realSubTimepointVarDiff = squeeze(realSubTimepointVar(:,:,2) - realSubTimepointVar(:,:,1));
permSubTimepointVar = squeeze(mean(permSubStd,5));
permSubTimepointVarDiff = squeeze(permSubTimepointVar(:,:,2,:) - permSubTimepointVar(:,:,1,:));
for iSub=1:length(goodSubs)
    for iRoi= ROIs
        pVal_sub_timepoint(iSub,iRoi) =  sum(permSubTimepointVarDiff(iSub,iRoi,:) >= realSubTimepointVarDiff(iSub,iRoi))/nperms;
    end
end


hypoth = pVal_sub<0.05;
hypoth_ph = pVal_fftPh_sub<0.05;
pVal_sub(:,ROIs)
pVal_fftPh_sub(:,ROIs)

meanRealAmpStdDiff = squeeze(mean(realSubAmpStdDiff));%mean over subjects
meanRealPhStdDiff = squeeze(mean(realSubPhStdDiff));%mean over subjects
meanRealOtherStdDiff = squeeze(mean(realSubOtherStdDiff));%mean over subjects
meanRealBaselineStdDiff = squeeze(mean(realSubBaselineStdDiff));%mean over subjects
meanRealFftAmpDiff = squeeze(mean(realSubMeanFftAmpDiff));%mean over subjects
meanRealFftPhDiff = squeeze(mean(realSubMeanFftPhDiff));%mean over subjects
meanRealStdVarDiff = squeeze(mean(realSubStdVarDiff));%mean over subjects

meanRealMeanStdDiff = squeeze(mean(realSubMeanStdDiff));%mean trials amplitude, mean over subjects
meanRealMedianStdDiff = squeeze(mean(realSubMedianStdDiff));%median trials amplitude, mean over subjects


meanPermAmpStdDiff = squeeze(mean(permSubAmpStdDiff));%mean over subjects
meanPermPhStdDiff = squeeze(mean(permSubPhStdDiff));%mean over subjects
meanPermOtherStdDiff = squeeze(mean(permSubOtherStdDiff));%mean over subjects
meanPermBaselineStdDiff = squeeze(mean(permSubBaselineStdDiff));%mean over subjects
meanPermMeanAmpDiff = squeeze(mean(permSubMeanAmpDiff));%mean over subjects
meanPermMeanPhDiff = squeeze(mean(permSubMeanPhDiff));%mean over subjects
meanPermMeanStdDiff = squeeze(mean(permSubMeanStdDiff));%mean over subjects
meanPermMedianStdDiff= squeeze(mean(permSubMedianStdDiff));%mean over subjects

meanPermStdVarDiff = squeeze(mean(permSubStdVarDiff));%mean over subjects


for iRoi= ROIs
    %RFX
    meanPermDiff(iRoi,:) = mean(permSubDiff(:,iRoi,:));%average over subjects
    meanRealDiff(iRoi) = mean(realSubDiff(:,iRoi));%average over subjects
    pVal_rfx(iRoi) = sum(meanPermDiff(iRoi,:) >= meanRealDiff(iRoi))/nperms;
    
    %FFT Phase RFX
    meanPermPhDiff(iRoi,:) = mean(permSubFftPhDiff(:,iRoi,:));%average over subjects
    meanRealPhDiff(iRoi) = mean(realSubFftPhDiff(:,iRoi));%average over subjects
    pVal_ph_rfx(iRoi) = sum(meanPermPhDiff(iRoi,:) >= meanRealPhDiff(iRoi))/nperms;
    
    meanPermAmpDiff(iRoi,:) = mean(permSubFftAmpDiff(:,iRoi,:));%average over subjects
    meanRealAmpDiff(iRoi) = mean(realSubFftAmpDiff(:,iRoi));%average over subjects
    pVal_amp_rfx(iRoi) = sum(meanPermAmpDiff(iRoi,:) >= meanRealAmpDiff(iRoi))/nperms;
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
    pVal_ffx(iRoi) = sum(meanPermDiffStd(iRoi,:) >= meanRealDiffStd(iRoi))/nperms;
    
    %variability
    permStd(iRoi,:,:,:) = mean(permSubStd(:,iRoi,:,:,:));%mean over subjects. permStd(iRoi,iRoi,rwd,p)
    for t=1:trialLength
       pVal_var_timecourse(iRoi,t) = sum( (permStd(iRoi,2,:,t)-permStd(iRoi,1,:,t)) >= (meanTimepointStd(iRoi,2,t)-meanTimepointStd(iRoi,1,t)))/nperms;
    end
    pVal_var(iRoi) = sum( mean(permStd(iRoi,2,:,:)-permStd(iRoi,1,:,:),4) >= mean(meanTimepointStd(iRoi,2,:)-meanTimepointStd(iRoi,1,:),3) )/nperms;
    
    %FFT variability
    pVal_ampVar(iRoi) = sum(meanPermAmpStdDiff(iRoi,:) >= meanRealAmpStdDiff(iRoi))/nperms;
    pVal_phVar(iRoi) = sum(meanPermPhStdDiff(iRoi,:) >= meanRealPhStdDiff(iRoi))/nperms;
    pVal_otherVar(iRoi) = sum(meanPermOtherStdDiff(iRoi,:) >= meanRealOtherStdDiff(iRoi))/nperms;
    pVal_baselineVar(iRoi) = sum(meanPermBaselineStdDiff(iRoi,:) >= meanRealBaselineStdDiff(iRoi))/nperms;

    pVal_phMean(iRoi) = sum(meanPermMeanPhDiff(iRoi,:) >= meanRealFftPhDiff(iRoi))/nperms;
    pVal_ampMean(iRoi) = sum(meanPermMeanAmpDiff(iRoi,:) >= meanRealFftAmpDiff(iRoi))/nperms;
    pVal_stdMean(iRoi) = sum(meanPermMeanStdDiff(iRoi,:) >= meanRealMeanStdDiff(iRoi))/nperms;
    pVal_stdMedian(iRoi) = sum(meanPermMedianStdDiff(iRoi,:) >= meanRealMedianStdDiff(iRoi))/nperms;
    
    pVal_stdVar(iRoi) = sum(meanPermStdVarDiff(iRoi,:) >= meanRealStdVarDiff(iRoi))/nperms;
    
    
end
pVal_rfx
pVal_ph_rfx
pVal_amp_rfx
% pVal_fft_rfx
pVal_var
% pVal_var_timecourse*trialLength

pVal_ampVar
pVal_phVar
pVal_otherVar
pVal_baselineVar
pVal_ampMean %mean trial amplitude measured by fft (per trial), low minus high
pVal_stdMean %mean trial amplitude measured by std (per trial), low minus high
pVal_stdMedian
pVal_phMean

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
realStdRT = std(realSubRT);
realRTdiff = realMeanRT(2) - realMeanRT(1);
pVal_RT = sum(permRTdiff >= realRTdiff)/nperms;



%%
histbins=1000;
rows=length(roiNames);
cols = length(goodSubs);


%%
% globalMean{iSub,rwd}\roiTC{iSub,iRoi,rwd}.tSeries';



clear smallSubAmp
for iRoi = 1:length(roiNames)
    for s=1:length(goodSubs)
        iSub = goodSubs(s);
        for rwd=1:2
            smallSubAmp(s,iRoi,rwd) = std(squeeze(subResponse(iSub,iRoi,rwd,:)));
        end
    end
end



smallSubDiff = smallSubAmp(:,:,1) - smallSubAmp(:,:,2);
subStdDiff = mean(subTimepointStd(:,:,2,:),4) - mean(subTimepointStd(:,:,1,:),4);


pVal_RTvar
realSubRTvarDiff = realSubRTvar(:,2) - realSubRTvar(:,1);

pVal_RT
realSubRTdiff = realSubRT(:,2)-realSubRT(:,1);

pVal_stdVar %variability of amplitude measured by std





%%
iRoi=2;
stdOfMeanTrial = mean(std(realSubMeanTC(:,iRoi,:,:),0,4));
meanOfTrialStd = mean(subMeanStd(:,iRoi,:));
meanOfTrialsFFTamp = mean(subFftMeanAmp(:,iRoi,:));
medianOfTrialsStd = subMedianStd(:,iRoi,:);

 temp = std(realSubMeanTC(:,iRoi,:,:),0,4);
 stdMeanTrialDiff = squeeze(temp(:,:,1) - temp(:,:,2));

 medianTrialStdDiff = squeeze(subMedianStd(:,iRoi,1) - subMedianStd(:,iRoi,2));
 
 %% MEAN SPECTRUM
iRoi=2;
 i=i+1; figure; clf
 TR=1.5;
 dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;
 for rwd=1:2
     runFFT = squeeze(mean(runMeanFFT(:,iRoi,rwd,:)));
     runFFTstd = squeeze(std(runMeanFFT(:,iRoi,rwd,:)));
     nframes = length(runFFT);
        %from plotMeanFourierAmp.m
        runFFT = runFFT / (nframes/2);
        frequencies = [0:nframes-1]/(nframes*TR);
        runFFT = runFFT(1:floor(nframes/2));
        frequencies = frequencies(1:floor(nframes/2));
        runFFTstd = runFFTstd / (nframes/2);
        runFFTstd = runFFTstd(1:floor(nframes/2));
        
%         plot(frequencies, runFFT,'.-','Color',plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);
%         
%         hold on
        
        dsErrorsurface(frequencies, runFFT, runFFTstd./sqrt(size(runMeanFFT,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(frequencies, runFFT,'.-', 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',markersize);

        
 end
 
 %%
 %% Timepoint variability as function of time
subTimepointStdDiff = squeeze(subTimepointStd(:,:,1,:) - subTimepointStd(:,:,2,:));
 i=i+1; figure; clf
 cols=5;
 rows=ceil(length(goodSubs)/cols);
 for iSub=1:length(goodSubs)
     subplot(rows,cols,iSub)
     for rwd=1:2
         plot(squeeze(subTimepointStd(iSub,iRoi,rwd,:)),'color', plotColors{rwd},'linewidth',2);
         hold on
     end
     plot(squeeze(subTimepointStdDiff(iSub,iRoi,:)),'k','linewidth',2);
 end
 groupMeanVar = mean(subTimepointStd(:,iRoi,:,:));
 subplot(rows,cols,rows*cols)
 for rwd=1:2
     plot(squeeze(groupMeanVar(:,:,rwd,:)),'color', plotColors{rwd},'linewidth',3);
     hold on
 end
 plot(squeeze(groupMeanVar(:,:,2,:)-groupMeanVar(:,:,1,:)),'k','linewidth',2);
  title('mean timepoint variability');
  

 
 for iRoi = 1:length(roiNames)
     for iSub=1:length(goodSubs)
         meanResp = squeeze(mean(subResponse(goodSubs(iSub),iRoi,:,:)));
         respDeriv(iSub,iRoi,1:9) = abs(diff(meanResp));
         respDeriv(iSub,iRoi,10) = abs(meanResp(end) - meanResp(1));
         varDiffRespDerivCorr(iSub,iRoi) = corr(squeeze(subTimepointStdDiff(iSub,iRoi,:)),squeeze(respDeriv(iSub,iRoi,:)));
     end
 end
 groupDeriv(1:9) = abs(diff(groupMeanVar(:,:,rwd,:)));
  groupDeriv(10) = abs(groupMeanVar(:,:,rwd,end) - groupMeanVar(:,:,rwd,1));
 corr(squeeze(groupMeanVar(:,:,2,:)-groupMeanVar(:,:,1,:)), squeeze(groupDeriv)')

%% correlate trial-by-trial RT with response amplitude/timing
%correlate RT with BOLD timing


for iRoi= ROIs
    for rwd=1:2
        allRTvec{iRoi,rwd}=[];
        allPhVec{iRoi,rwd}=[];
        allStdVec{iRoi,rwd}=[];
        allAmpVec{iRoi,rwd}=[];
        for iSub=1:length(goodSubs)
            [rtBoldPhCorr(iSub,iRoi,rwd) rtBoldPhCorrPval(iSub,iRoi,rwd)]=  corr(trialRoiFftPhVec{iSub,iRoi,rwd}', trialRTvec{iSub,rwd});
            [rtBoldPhShiftCorr(iSub,iRoi,rwd) rtBoldPhShiftCorrPval(iSub,iRoi,rwd)]=  corr(circshift(trialRoiFftPhVec{iSub,iRoi,rwd}',-1), trialRTvec{iSub,rwd});
            
            [rtBoldAmpCorr(iSub,iRoi,rwd) rtBoldAmpCorrPval(iSub,iRoi,rwd)]=  corr(trialRoiFftAmpVec{iSub,iRoi,rwd}', trialRTvec{iSub,rwd});
            [rtBoldAmpShiftCorr(iSub,iRoi,rwd) rtBoldAmpShiftCorrPval(iSub,iRoi,rwd)]=  corr(circshift(trialRoiFftAmpVec{iSub,iRoi,rwd}',-1), trialRTvec{iSub,rwd});
            
            [rtBoldStdCorr(iSub,iRoi,rwd) rtBoldStdCorrPval(iSub,iRoi,rwd)]=  corr(trialRoiAmpVec{iSub,iRoi,rwd}', trialRTvec{iSub,rwd});
            [rtBoldStdShiftCorr(iSub,iRoi,rwd) rtBoldStdShiftCorrPval(iSub,iRoi,rwd)]=  corr(circshift(trialRoiAmpVec{iSub,iRoi,rwd}',-1), trialRTvec{iSub,rwd});
            
            
            subTrialMat = [trialRoiFftPhVec{iSub,iRoi,rwd}', trialRoiFftAmpVec{iSub,iRoi,rwd}',trialRoiAmpVec{iSub,iRoi,rwd}',-trialRTvec{iSub,rwd}];
            [subTrialCorr(iSub,iRoi,rwd,:,:) subTrialCorrPval(iSub,iRoi,rwd,:,:)] = corr(subTrialMat);
            
            %pool across subjects
%             allRTvec{iRoi,rwd} = [allRTvec{iRoi,rwd}; zscore(trialRTvec{iSub,rwd})]; 
%             allPhVec{iRoi,rwd}= [allPhVec{iRoi,rwd}; zscore(trialRoiFftPhVec{iSub,iRoi,rwd})'];
%             allStdVec{iRoi,rwd}=[allStdVec{iRoi,rwd}; zscore(trialRoiAmpVec{iSub,iRoi,rwd})'];
%             allAmpVec{iRoi,rwd}=[allAmpVec{iRoi,rwd}; zscore(trialRoiFftAmpVec{iSub,iRoi,rwd})'];
            
            allRTvec{iRoi,rwd} = [allRTvec{iRoi,rwd}; trialRTvec{iSub,rwd}]; 
            allPhVec{iRoi,rwd}= [allPhVec{iRoi,rwd}; trialRoiFftPhVec{iSub,iRoi,rwd}'];
            allStdVec{iRoi,rwd}=[allStdVec{iRoi,rwd}; trialRoiAmpVec{iSub,iRoi,rwd}'];
            allAmpVec{iRoi,rwd}=[allAmpVec{iRoi,rwd}; trialRoiFftAmpVec{iSub,iRoi,rwd}'];
        end
        rwdTrialMat{iRoi, rwd} = [allPhVec{iRoi,rwd}'; allAmpVec{iRoi,rwd}'; allStdVec{iRoi,rwd}'; -allRTvec{iRoi,rwd}'];
        [roiTrialCorr(iRoi,rwd,:,:) roiTrialCorrPval(iRoi,rwd,:,:)] = corr(rwdTrialMat{iRoi, rwd}');
    end
    trialMat{iRoi} = [rwdTrialMat{iRoi, 1} rwdTrialMat{iRoi, 2}]; 
    [roiTrialCorr(iRoi,3,:,:) roiTrialCorrPval(iRoi,3,:,:)] = corr(trialMat{iRoi}');
end
trialLabels = {'FFT ph','FFT amp','STD amp','-RT'};

%% pooled trial correlations
i=i+1; figure; clf
iRoi=2;
rows=3;
cols=3;
for rwd=1:2
    subplot(rows,cols,rwd)
    imagesc(squeeze( roiTrialCorr(iRoi,rwd,:,:)));
    yticks(1:length(trialLabels));
    yticklabels(trialLabels);
    xticks(1:length(trialLabels));
    xticklabels(trialLabels);
    title(['pooled trial-by-trial correlations, ' rwdString{rwd}]);
end
subplot(rows,cols,3)
imagesc(squeeze( roiTrialCorr(iRoi,3,:,:)));


%permutation test of pooled correlations
nfactors = size(rwdTrialMat{iRoi,rwd},1);
ntotaltrials = size(rwdTrialMat{iRoi,rwd},2);
clear permTrialMat
for iRoi=ROIs
    for rwd=1:2
        for p=1:nperms
            clear permRwdTrialMat
            for ifactor=1:nfactors
                permRwdTrialMat(ifactor,:) = rwdTrialMat{iRoi,rwd}(ifactor,randperm(ntotaltrials));
            end
            permTrialCorr(iRoi,rwd,p,:,:) = corr(permRwdTrialMat');
        end
    end
    for p=1:nperms
        clear permTrialMat
        for ifactor=1:nfactors
            permTrialMat(ifactor,:) = trialMat{iRoi}(ifactor,randperm(ntotaltrials));
        end
        permTrialCorr(iRoi,3,p,:,:) = corr(permTrialMat');
    end
end

%p-values for trial-by-trial correlations
for iRoi=ROIs
    for rwd=1:3
        for ifactor=1:nfactors
            for jfactor=1:nfactors
                pvalTrialCorr(iRoi,rwd,ifactor,jfactor) = sum(squeeze(permTrialCorr(iRoi,rwd,:,ifactor,jfactor))>=roiTrialCorr(iRoi,rwd,ifactor,jfactor))./nperms;
            end
        end
    end
end

for rwd=1:3
    subplot(rows,cols,cols+rwd)
    imagesc(squeeze( pvalTrialCorr(2,rwd,:,:)));
end


for rwd=1:3
    subplot(rows,cols,2*cols+rwd)
    imagesc(squeeze( roiTrialCorrPval(2,rwd,:,:)));
end

%% mean over single subject correlations
for iRoi=ROIs
    for rwd=1:2
        meanSubTrialCorr(iRoi,rwd,:,:) = squeeze(mean(subTrialCorr(:,iRoi,rwd,:,:)));
    end
end
i=i+1; figure; clf
 rows=1; 
 cols=2;
 for rwd=1:2
     subplot(rows,cols,rwd)
     imagesc(squeeze(meanSubTrialCorr(2,rwd,:,:)));
     yticks(1:length(trialLabels));
yticklabels(trialLabels);
xticks(1:length(trialLabels));
xticklabels(trialLabels);
title(['subject trial-by-trial correlations, ' rwdString{rwd}]);
 end
 
%% CORRELATION MATRIX
iRoi=2;
% analysesMat = [smallSubDiff(:,iRoi),realSubMeanFftPhDiff(:,iRoi),subStdDiff(:,iRoi),realSubPhStdDiff(:,iRoi),realSubAmpStdDiff(:,iRoi),realSubStdVarDiff(:,iRoi),realSubRTvarDiff,realSubRTdiff];
% lbls = {'std amp','mean phase', 'timepoint var','FFT ph var','FFT amp var','std amp var','RT var','mean RT' };
% analysesMat = [smallSubDiff(:,iRoi),realSubMeanFftPhDiff(:,iRoi),realSubPhStdDiff(:,iRoi),realSubAmpStdDiff(:,iRoi),realSubStdVarDiff(:,iRoi),subStdDiff(:,iRoi),realSubRTvarDiff,realSubRTdiff];
% lbls = {'std amp','mean phase','FFT ph var','FFT amp var','std amp var', 'timepoint var','RT var','mean RT' };

analysesMat = [realSubMeanFftPhDiff(:,iRoi),realSubPhStdDiff(:,iRoi),smallSubDiff(:,iRoi),realSubAmpStdDiff(:,iRoi),realSubStdVarDiff(:,iRoi),subStdDiff(:,iRoi),realSubRTvarDiff,realSubRTdiff];
lbls = {'mean phase','FFT ph var','std amp','FFT amp var','std amp var', 'timepoint var','RT var','mean RT' };


[rmat pmat] = corr(analysesMat);
for ianalysis=1:size(pmat,1)
%     rmat(ianalysis,ianalysis) = 1;
    pmat(ianalysis,ianalysis) = 0;
end

 i=i+1; figure; clf
 rows=1; 
 cols=2;
subplot(rows,cols,1)
imagesc(rmat)

yticks(1:length(lbls));
yticklabels(lbls);

 subplot(rows,cols,2)
 imagesc(pmat<0.05);
yticks(1:length(lbls));
yticklabels(lbls);

i=i+1; figure; clf
pd = pdist(analysesMat','correlation');
Y = cmdscale(pd);
% Y = mdscale(pd,2);
scatter(Y(:,1),Y(:,2),70,1:size(Y,1))
colormap jet
%%
mean(subPropCorrect)
std(subPropCorrect)
mean(realSubRT)
std(realSubRT)