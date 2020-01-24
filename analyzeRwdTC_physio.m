close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

physioRandSeed = rng;
% load(fullfile(dataFolder,'randSeed.mat'),'randSeed');

nperms=1000;
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
if regressGlobalMean
    globalMeanString = '_globalRegressed';
end
ConcatProjStr = '';
if ConcatProj
    ConcatProjStr = 'ConcatProj';
end
load([dataFolder 'rwdTC_physio' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'globalMean','regressBetasGlobal','runRwd',...
    'subRoiRuns','runMeanFFT',...
    'allVoxTrialResponse','allVoxTaskPhase','allVoxTaskAmp','allVoxTaskCo',...
    'voxTrials','voxGoodTrials','meanVoxTrial',...
    'maxRT',...
    'ecgselect','ecgSampleRate','ecgTrial','ecgRunLength','ecgInterpMethod',...
    'ecg','ecgPulseRate','interpPulseRate',...
    'respselect','resp',...
    'rwdPulseTC','rwdRvTC',...
    'designMatPulse','designMatRespPulse','designMatResp','deconvLength',...
    'allGoodTrials');

plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);

goodSubs = [1:length(subFolders)]; 

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



        trRV{iSub,rwd} = reshape(rwdRvTC{iSub,rwd}, ecgTrial,[]);%this is only good trials!
        subMeanRV(iSub,rwd,:) = nanmean(trRV{iSub,rwd},2);
        subRespStd(iSub,rwd) = std(subMeanRV(iSub,rwd,:));%std amplitude of mean
        subRespVar(iSub,rwd,:) = nanstd(trRV{iSub,rwd},0,2);%timepoint variability
        respStd = nanstd(trRV{iSub,rwd});%std amp per trial
        subRespStdVar(iSub,rwd) = std(respStd);
        f = fft(trRV{iSub,rwd});
        subRespPhVar(iSub,rwd) = nanstd(angle(f(2,:)));
        subRespAmpVar(iSub,rwd) = nanstd(abs(f(2,:)));
        
        
        trPulse{iSub,rwd} = reshape(rwdPulseTC{iSub,rwd}, ecgTrial,[]);%this is only good trials!
        subMeanPulse(iSub,rwd,:) = nanmean(trPulse{iSub,rwd},2);
        subPulseStd(iSub,rwd) = std(subMeanPulse(iSub,rwd,:));%std amplitude of mean
        subPulseVar(iSub,rwd,:) = nanstd(trPulse{iSub,rwd},0,2);%timepoint variability
        pulseStd = nanstd(trPulse{iSub,rwd});%std amp per trial
        subPulseStdVar(iSub,rwd) = std(pulseStd);
        f = fft(trPulse{iSub,rwd});
        subPulsePhVar(iSub,rwd) = nanstd(angle(f(2,:)));
        subPulseAmpVar(iSub,rwd) = nanstd(abs(f(2,:)));
        
            %REGRESS OUT PHYSIO
            %roiMeanTseries{iSub,iRoi,rwd} includes all trials
           pulseKernel(iSub,iRoi,rwd,:) = designMatPulse{iSub,rwd}'\subTrialResponse{iSub,iRoi,rwd}(:);%this is only good trials!
           pulseResidualTC{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:)' - squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{iSub,rwd};

           respKernel(iSub,iRoi,rwd,:) = designMatResp{iSub,rwd}'\subTrialResponse{iSub,iRoi,rwd}(:);
           respResidualTC{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:)' - squeeze(respKernel(iSub,iRoi,rwd,:))'*designMatResp{iSub,rwd};

           physioKernel(iSub,iRoi,rwd,:) = designMatRespPulse{iSub,rwd}'\subTrialResponse{iSub,iRoi,rwd}(:);
           physioResidualTC{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:)' - squeeze(physioKernel(iSub,iRoi,rwd,:))'*designMatRespPulse{iSub,rwd};

           
%            subTrialResponse{iSub,iRoi,rwd} = reshape(roiMeanTseries{iSub,iRoi,rwd}(:), trialLength, length(roiMeanTseries{iSub,iRoi,rwd}(:))/trialLength);

            %KEEP ONLY THE RESIDUAL
%            subTrialResponse{iSub,iRoi,rwd} = reshape(pulseResidualTC{iSub,iRoi,rwd}(:), trialLength, length(pulseResidualTC{iSub,iRoi,rwd}(:))/trialLength);
%            subTrialResponse{iSub,iRoi,rwd} = reshape(respResidualTC{iSub,iRoi,rwd}(:), trialLength, length(respResidualTC{iSub,iRoi,rwd}(:))/trialLength);
            subTrialResponse{iSub,iRoi,rwd} = reshape(physioResidualTC{iSub,iRoi,rwd}(:), trialLength, length(physioResidualTC{iSub,iRoi,rwd}(:))/trialLength);

           
           reshapedTrials = subTrialResponse{goodSubs(iSub),iRoi,rwd};
           temp = fft(reshapedTrials);
           roiFftAmp = abs(temp(2,:));
           roiFftBaseline = abs(temp(1,:));

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
    
    subRV = [trRV{iSub,1} trRV{iSub,2}];%750 timepoints, trials
%     subRespVar = nanstd(subRV);%std amp per trial
    subPulse = [trPulse{iSub,1} trPulse{iSub,2}];
%     subPulseVar = nanstd(subPulse);%std amp per trial
    
    %pulse & respiration deconvolution design matrices
    temp = [designMatPulse{iSub,1} designMatPulse{iSub,2}];
    designMatPulseTrials = reshape(temp,deconvLength,trialLength,[]);%10 deconvolution points, 10 trial timepoints, X trials
    temp = [designMatResp{iSub,1} designMatResp{iSub,2}];
    designMatRespTrials = reshape(temp,deconvLength,trialLength,[]);%10 deconvolution points, 10 trial timepoints, X trials
    temp = [designMatRespPulse{iSub,1} designMatRespPulse{iSub,2}];
    designMatPhysioTrials = reshape(temp,deconvLength,trialLength,[]);%10 deconvolution points, 10 trial timepoints, X trials
 
        
    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));
        

        

        for rwd=1:2
            %mean pulse + respiration
            permRV = subRV(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));
            permSubRespStd(iSub,rwd,p) = std(nanmean(permRV,2));%std amp of avg
            permSubRespStdVar(iSub,rwd,p) = std(nanstd(permRV));%std amp variability
%             permSubRespStdVar(iSub,rwd,p) = std(subRespVar(randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1)));%std amp variability
            permSubMeanRV(iSub,rwd,p,:) = nanmean(permRV,2);%mean  per permutation
            permSubRespVar(iSub,rwd,p,:) = nanstd(permRV,0,2);%timepoint variability
            f = fft(permRV);
            permSubRespPhVar(iSub,rwd,p) = nanstd(angle(f(2,:)));
            permSubRespAmpVar(iSub,rwd,p) = nanstd(abs(f(2,:)));
            
            
            
            permPulse = subPulse(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));
            permSubPulseStd(iSub,rwd,p) = std(nanmean(permPulse,2));%std amp of avg
            permSubPulseStdVar(iSub,rwd,p) = std(nanstd(permPulse));%std amp variability
%             permSubPulseStdVar(iSub,rwd,p) = std(subPulseVar(randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1)));%std amp variability
            permSubMeanPulse(iSub,rwd,p,:) = nanmean(permPulse,2);
            permSubPulseVar(iSub,rwd,p,:) = nanstd(permRV,0,2);%timepoint variability
            f = fft(permPulse);
            permSubPulsePhVar(iSub,rwd,p) = nanstd(angle(f(2,:)));
            permSubPulseAmpVar(iSub,rwd,p) = nanstd(abs(f(2,:)));

            
            %physio deconvolution kernels
            temp = designMatPulseTrials(:,:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));
            permDesignMatPulse{iSub,rwd} = reshape(temp,deconvLength,[]);
            temp = designMatRespTrials(:,:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));
            permDesignMatResp{iSub,rwd} = reshape(temp,deconvLength,[]);
            temp = designMatPhysioTrials(:,:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));
            permDesignMatPhysio{iSub,rwd} = reshape(temp,deconvLength,[]);

        end
        
        
        
        for iRoi= ROIs%length(roiNames)
            for rwd=1:2
                permTrials = roiTrials{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));% 10 timepoints X #trials
                permSubMeanTC(iSub,iRoi,rwd,p,:) = mean(permTrials,2);%average timecourse over trials
                permSubResp(iSub,iRoi,rwd,p) = std(squeeze(permSubMeanTC(iSub,iRoi,rwd,p,:)));%amplitude of task response
                
                %latency
                temp=fft(permSubMeanTC(iSub,iRoi,rwd,p,:));
                permSubMeanFftAmp(iSub,iRoi,rwd,p) = abs(temp(2));
                permSubMeanFftPh(iSub,iRoi,rwd,p) = angle(temp(2));

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
                
                %physio
                permPulseKernel(iSub,iRoi,rwd,p,:) = permDesignMatPulse{iSub,rwd}'\permTrials(:);
                permRespKernel(iSub,iRoi,rwd,p,:) = permDesignMatResp{iSub,rwd}'\permTrials(:);
                permPhysioKernel(iSub,iRoi,rwd,p,:) = permDesignMatPhysio{iSub,rwd}'\permTrials(:);
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


% permSubRespStd = nanstd(permSubMeanRV,0,4);

%permSubPulseStd permSubMeanPulse permSubPulseVar permSubPulsePhVar permSubPulseAmpVar
%subPulseStd subPulseVar subPulsePhVar subPulseAmpVar


permSubRespStdDiff = squeeze(permSubRespStd(:,1,:) - permSubRespStd(:,2,:));%iSub,perm
permSubRespStdVarDiff = squeeze(permSubRespStdVar(:,2,:) - permSubRespStdVar(:,1,:));
permSubRespPhVarDiff =  squeeze(permSubRespPhVar(:,2,:) - permSubRespPhVar(:,1,:));
permSubRespAmpVarDiff =  squeeze(permSubRespAmpVar(:,2,:) - permSubRespAmpVar(:,1,:));
permSubRespVarMean = squeeze(mean(permSubRespVar,4));
permSubRespVarDiff = squeeze(permSubRespVarMean(:,2,:) - permSubRespVarMean(:,1,:));

permSubPulseStdDiff = squeeze(permSubPulseStd(:,1,:) - permSubPulseStd(:,2,:));
permSubPulseStdVarDiff = squeeze(permSubPulseStdVar(:,2,:) - permSubPulseStdVar(:,1,:));
permSubPulsePhVarDiff =  squeeze(permSubPulsePhVar(:,2,:) - permSubPulsePhVar(:,1,:));
permSubPulseAmpVarDiff =  squeeze(permSubPulseAmpVar(:,2,:) - permSubPulseAmpVar(:,1,:));
permSubPulseVarMean = squeeze(mean(permSubPulseVar,4));
permSubPulseVarDiff = squeeze(permSubPulseVarMean(:,2,:) - permSubPulseVarMean(:,1,:));

%PULSE P VALUES
realSubPulseStdDiff = subPulseStd(:,1) - subPulseStd(:,2);
pval_pulseStd = sum(mean(permSubPulseStdDiff)>=mean(realSubPulseStdDiff))/nperms;

realSubPulseStdVarDiff = subPulseStdVar(:,2) - subPulseStdVar(:,1);
pval_pulseStdVar = sum(mean(permSubPulseStdVarDiff)>=mean(realSubPulseStdVarDiff))/nperms;

realSubPulsePhVarDiff = subPulsePhVar(:,2) - subPulsePhVar(:,1);
pval_pulsePhVar = sum(mean(permSubPulsePhVarDiff)>=mean(realSubPulsePhVarDiff))/nperms;

realSubPulseAmpVarDiff = subPulseAmpVar(:,2) - subPulseAmpVar(:,1);
pval_pulseAmpVar = sum(mean(permSubPulseAmpVarDiff)>=mean(realSubPulseAmpVarDiff))/nperms;

subPulseVarMean = mean(subPulseVar,3);%mean timepoint variability
realSubPulseVarMeanDiff = subPulseVarMean(:,2) - subPulseVarMean(:,1);
pval_pulseVarMean = sum(mean(permSubPulseVarDiff)>=mean(realSubPulseVarMeanDiff))/nperms;

%RESPIRATION P VALUES
realSubRespStdDiff = subRespStd(:,1) - subRespStd(:,2);
pval_respStd = sum(mean(permSubRespStdDiff)>=mean(realSubRespStdDiff))/nperms;

realSubRespStdVarDiff = subRespStdVar(:,2) - subRespStdVar(:,1);
pval_respStdVar = sum(mean(permSubRespStdVarDiff)>=mean(realSubRespStdVarDiff))/nperms;

realSubRespPhVarDiff = subRespPhVar(:,2) - subRespPhVar(:,1);
pval_respPhVar = sum(mean(permSubRespPhVarDiff)>=mean(realSubRespPhVarDiff))/nperms;

realSubRespAmpVarDiff = subRespAmpVar(:,2) - subRespAmpVar(:,1);
pval_respAmpVar = sum(mean(permSubRespAmpVarDiff)>=mean(realSubRespAmpVarDiff))/nperms;

subRespVarMean = mean(subRespVar,3);%mean timepoint variability
realSubRespVarMeanDiff = subRespVarMean(:,2) - subRespVarMean(:,1);
pval_respVarMean = sum(mean(permSubRespVarDiff)>=mean(realSubRespVarMeanDiff))/nperms;

keyboard



realSubPulseKernelStd = std(pulseKernel,0,4);
realSubPulseKernelStdDiff = realSubPulseKernelStd(:,:,1) - realSubPulseKernelStd(:,:,2);
meanRealPulseKernelStdDiff = squeeze(mean(realSubPulseKernelStdDiff));
permSubPulseKernelStd = std(permPulseKernel,0,5);%permPulseKernel(iSub,iRoi,rwd,p,:)
permSubPulseKernelStdDiff = permSubPulseKernelStd(:,:,1,:) - permSubPulseKernelStd(:,:,2,:);
meanPermPulseKernelStdDiff = squeeze(mean(permSubPulseKernelStdDiff));
pval_pulseKernel_std = sum(meanPermPulseKernelStdDiff'>=meanRealPulseKernelStdDiff)./nperms;

realSubRespKernelStd = std(respKernel,0,4);
realSubRespKernelStdDiff = realSubRespKernelStd(:,:,1) - realSubRespKernelStd(:,:,2);
meanRealRespKernelStdDiff = squeeze(mean(realSubRespKernelStdDiff));
permSubRespKernelStd = std(permRespKernel,0,5);%permPulseKernel(iSub,iRoi,rwd,p,:)
permSubRespKernelStdDiff = permSubRespKernelStd(:,:,1,:) - permSubRespKernelStd(:,:,2,:);
meanPermRespKernelStdDiff = squeeze(mean(permSubRespKernelStdDiff));
pval_respKernel_std = sum(meanPermRespKernelStdDiff'>=meanRealRespKernelStdDiff)./nperms;

% realSubPulseStd = std(pulseKernel,0,4);
% realSubPulseStdDiff = realSubPulseStd(:,:,1) - realSubPulseStd(:,:,2);
% meanRealPulseStdDiff = squeeze(mean(realSubPulseStdDiff));
% permSubPulseStd = std(permPulseKernel,0,5);%permPulseKernel(iSub,iRoi,rwd,p,:)
% permSubPulseStdDiff = permSubPulseStd(:,:,1,:) - permSubPulseStd(:,:,2,:);
% meanPermPulseStdDiff = squeeze(mean(permSubPulseStdDiff));
% pval_pulse_std = sum(meanPermPulseStdDiff'>=meanRealPulseStdDiff)./nperms;
% 
% realSubRespStd = std(respKernel,0,4);
% realSubRespStdDiff = realSubRespStd(:,:,1) - realSubRespStd(:,:,2);
% meanRealRespStdDiff = squeeze(mean(realSubRespStdDiff));
% permSubRespStd = std(permRespKernel,0,5);%permPulseKernel(iSub,iRoi,rwd,p,:)
% permSubRespStdDiff = permSubRespStd(:,:,1,:) - permSubRespStd(:,:,2,:);
% meanPermRespStdDiff = squeeze(mean(permSubRespStdDiff));
% pval_resp_std = sum(meanPermRespStdDiff'>=meanRealRespStdDiff)./nperms;

% realSubPhysioStd = std(physioKernel,0,4);
% realSubPhysioStdDiff = realSubPhysioStd(:,:,1) - realSubPhysioStd(:,:,2);
% meanRealPhysioStdDiff = squeeze(mean(realSubPhysioStdDiff));
% permSubPhysioStd = std(permPhysioKernel,0,5);%permPulseKernel(iSub,iRoi,rwd,p,:)
% permSubPhysioStdDiff = permSubPhysioStd(:,:,1,:) - permSubPhysioStd(:,:,2,:);
% meanPermPhysioStdDiff = squeeze(mean(permSubPhysioStdDiff));
% pval_physio_std = sum(meanPermPhysioStdDiff'>meanRealPhysioStdDiff);


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


save(fullfile(dataFolder,'physioRandSeed.mat'),'physioRandSeed');


%%
globalMean{iSub,rwd}\roiTC{iSub,iRoi,rwd}.tSeries';



clear smallSubAmp
for iRoi = 1:length(roiNames)
    for s=1:length(goodSubs)
        iSub = goodSubs(s);
        for rwd=1:2
            smallSubAmp(s,iRoi,rwd) = std(squeeze(subResponse(iSub,iRoi,rwd,:)));
        end
    end
end


iRoi=2;

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
  


%% correlate trial-by-trial RT with response amplitude/timing
%correlate RT with BOLD timing

for iSub=1:length(goodSubs)
    for iRoi= ROIs
        for rwd=1:2
            
           [rtBoldPhCorr(iSub,iRoi,rwd) rtBoldPhCorrPval(iSub,iRoi,rwd)]=  corr(trialRoiFftPhVec{iSub,iRoi,rwd}', trialRTvec{iSub,rwd});
           [rtBoldPhShiftCorr(iSub,iRoi,rwd) rtBoldPhShiftCorrPval(iSub,iRoi,rwd)]=  corr(circshift(trialRoiFftPhVec{iSub,iRoi,rwd}',-1), trialRTvec{iSub,rwd});
            
           [rtBoldAmpCorr(iSub,iRoi,rwd) rtBoldAmpCorrPval(iSub,iRoi,rwd)]=  corr(trialRoiFftAmpVec{iSub,iRoi,rwd}', trialRTvec{iSub,rwd});
            [rtBoldAmpShiftCorr(iSub,iRoi,rwd) rtBoldAmpShiftCorrPval(iSub,iRoi,rwd)]=  corr(circshift(trialRoiFftAmpVec{iSub,iRoi,rwd}',-1), trialRTvec{iSub,rwd});

            [rtBoldStdCorr(iSub,iRoi,rwd) rtBoldStdCorrPval(iSub,iRoi,rwd)]=  corr(trialRoiAmpVec{iSub,iRoi,rwd}', trialRTvec{iSub,rwd});
            [rtBoldStdShiftCorr(iSub,iRoi,rwd) rtBoldStdShiftCorrPval(iSub,iRoi,rwd)]=  corr(circshift(trialRoiAmpVec{iSub,iRoi,rwd}',-1), trialRTvec{iSub,rwd});
            
            
            trialMat = [-trialRoiFftPhVec{iSub,iRoi,rwd}', trialRoiFftAmpVec{iSub,iRoi,rwd}',trialRoiAmpVec{iSub,iRoi,rwd}',trialRTvec{iSub,rwd}];
            [trialCorr(iSub,iRoi,rwd,:,:) trialCorrPval(iSub,iRoi,rwd,:,:)] = corr(trialMat);
        end
    end
end
trialLabels = {'FFT ph','FFT amp','STD amp','RT'};

for iRoi=ROIs
    for rwd=1:2
        meanTrialCorr(iRoi,rwd,:,:) = squeeze(mean(trialCorr(:,iRoi,rwd,:,:)));
    end
end
i=i+1; figure; clf
 rows=1; 
 cols=2;
 for rwd=1:2
     subplot(rows,cols,rwd)
     imagesc(squeeze(meanTrialCorr(2,rwd,:,:)));
     yticks(1:length(trialLabels));
yticklabels(trialLabels);
xticks(1:length(trialLabels));
xticklabels(trialLabels);
title(['trial-by-trial correlations, ' rwdString{rwd}]);
 end
 

%%
% mean(subPropCorrect)
% std(subPropCorrect)
% mean(realSubRT)
% std(realSubRT)

%%
pval_pulseKernel_std
pval_respKernel_std