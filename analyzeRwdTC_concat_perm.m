close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

nperms=10000;
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
clear subStd
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

% subject 4 = s0022
% subject 16 = outlier of RT and RT std
% goodSubs = [1:3 5:18];%excluding subject 22.     0.0830    0.0402
% goodSubs = [1:3 5:6 8:9 11 13:18];%excluding subject 22 and all double
% goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% subjects. 0.0457  0.0140,     0.0553    0.0214

for iSub = 1:length(goodSubs)%length(subdirs)
%     iSub = goodSubs(i)
     for iRoi=1:length(roiNames)
        for rwd=1:2
            reshapedTrials = reshape(subTrialResponse{goodSubs(iSub),iRoi,rwd},10,[]);
            if iSub==goodSubs(1)
                allTrials{iRoi,rwd} = reshapedTrials;
            else
                allTrials{iRoi,rwd} = [allTrials{iRoi,rwd} reshapedTrials];
            end
            
           %trial-by-trial variability per subject
           subStd(iSub,iRoi,rwd,:) = std(reshapedTrials,0,2);
           
           temp = fft(reshapedTrials);
           roiBinFftAmp = abs(temp(2,:));
           roiBinFftPh = angle(temp(2,:));
           roiBinFftOther = sum(abs(temp([1 3:end],:)));
           subBinFftAmpStd(iSub,iRoi,rwd) = std(roiBinFftAmp);
           subBinFftPhStd(iSub,iRoi,rwd) = circ_std(roiBinFftPh');
           subBinFftOtherStd(iSub,iRoi,rwd) = std(roiBinFftOther);
           roiBinFftBaseline = abs(temp(1,:));
           subBinFftBaselineStd(iSub,iRoi,rwd) = std(roiBinFftBaseline);
%            rehapedTrials = reshapedTrials - repmat(roiBinFftBaseline,10,1);
%            rehapedTrials = rehapedTrials./repmat(roiBinFftAmp,10,1);
           
        end
     end
end


for rwd=1:2
    for iRoi = 1:length(roiNames)
        meanResponse(iRoi,rwd,:) = mean(subResponse(goodSubs,iRoi,rwd,:));
        stdResponse(iRoi,rwd,:) = std(subResponse(goodSubs,iRoi,rwd,:));
%         plot(squeeze(meanResponse(rwd,iRoi,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
%         hold on
        meanStd(iRoi, rwd,:) = mean(subStd(:,iRoi,rwd,:));%mean trial-by-trial variability
        stdStd(iRoi, rwd,:) = std(subStd(:,iRoi,rwd,:));%std of trial-by-trial variability
    end
end


%% PERMUTATIONS
tic

for iSub = 1:length(goodSubs)%length(subdirs)
    for rwd=1:2

        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!

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
                           permFftBaseline = abs(temp(1,:));
                permSubBaselineStd(iSub,iRoi,rwd) = std(permFftBaseline);
            end
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

        
    end
end
hypoth = pVal_sub<0.05;
% pVal_sub(:,ROIs);

meanRealAmpStdDiff = squeeze(mean(realSubAmpStdDiff));%mean over subjects
meanRealPhStdDiff = squeeze(mean(realSubPhStdDiff));%mean over subjects
meanRealOtherStdDiff = squeeze(mean(realSubOtherStdDiff));%mean over subjects
meanPermAmpStdDiff = squeeze(mean(permSubAmpStdDiff));%mean over subjects
meanPermPhStdDiff = squeeze(mean(permSubPhStdDiff));%mean over subjects
meanPermOtherStdDiff = squeeze(mean(permSubOtherStdDiff));%mean over subjects


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
       pVal_var_timecourse(iRoi,t) = sum( (permStd(iRoi,2,:,t)-permStd(iRoi,1,:,t)) > (meanStd(iRoi,2,t)-meanStd(iRoi,1,t)))/nperms;
    end
    pVal_var(iRoi) = sum( mean(permStd(iRoi,2,:,:)-permStd(iRoi,1,:,:),4) > mean(meanStd(iRoi,2,:)-meanStd(iRoi,1,:),3) )/nperms;
    
    %FFT variability
    pVal_ampVar(iRoi) = sum(meanPermAmpStdDiff(iRoi,:) > meanRealAmpStdDiff(iRoi))/nperms;
    pVal_phVar(iRoi) = sum(meanPermPhStdDiff(iRoi,:) > meanRealPhStdDiff(iRoi))/nperms;
    pVal_otherVar(iRoi) = sum(meanPermOtherStdDiff(iRoi,:) > meanRealOtherStdDiff(iRoi))/nperms;
    
    
end
pVal_rfx
pVal_ffx
% pVal_fft_rfx
pVal_var
pVal_var_timecourse*trialLength

pVal_ampVar
pVal_phVar
pVal_otherVar

%% ANOVA
% varNames = {'baseline', 'std'};
varNames = {'ampVarH', 'ampVarL','phVarH','phVarL'};
varNames = {'Y1','Y2','Y3','Y4'};
% varNames = {'Y1','Y2'};
% tbl = table(realSubBaseline(:,1),realSubBaseline(:,2),realSubStd(:,1),realSubStd(:,2),'VariableNames',varNames);
iRoi=2;
ampVarH = subBinFftAmpStd(:,iRoi,1);
ampVarL = subBinFftAmpStd(:,iRoi,2);
phVarH = subBinFftPhStd(:,iRoi,1);
phVarL = subBinFftPhStd(:,iRoi,2);
tbl = array2table([ampVarH ampVarL phVarH phVarL], 'VariableNames',varNames);
reward = {'H','L','H','L'};
polarCoord = {'amplitude','amplitude','phase','phase'};

within = table(reward',polarCoord', 'VariableNames',{'reward','polarCoord'});

rm = fitrm(tbl, 'Y1-Y4~1', 'WithinDesign', within);
ranovatbl = ranova(rm, 'WithinModel', 'reward*polarCoord')


%%
histbins=1000;
rows=length(roiNames);
cols = length(goodSubs);
i=i+1;
figure(i)
clf
for iSub = 1:length(goodSubs)
    for iRoi= ROIs
        subplot(rows,cols,iSub + (iRoi-1)*cols)
        histogram(squeeze(permSubDiff(iSub,iRoi,:)),histbins); hold all
        vline(prctile(permSubDiff(iSub,iRoi,:),95),'k');
        vline(realSubDiff(iSub,iRoi),'r');
    end
end
        
        
%%
i=i+1;
figure(i)
clf
rows=length(roiNames);
cols = 2;
for iRoi= ROIs
    %RFX
    subplot(rows,cols,1+(iRoi-1)*cols)
        histogram(squeeze(meanPermDiff(iRoi,:)),histbins); hold all
        vline(prctile(meanPermDiff(iRoi,:),95),'k');
        vline(meanRealDiff(iRoi),'r');
        title(['random effect: ' roiNames{iRoi}]);
    %FFX
        subplot(rows,cols,2+(iRoi-1)*cols)
        histogram(squeeze(meanPermDiffStd(iRoi,:)),histbins); hold all
        vline(prctile(meanPermDiffStd(iRoi,:),95),'k');
        vline(meanRealDiffStd(iRoi),'r');
        title(['fixed effect: ' roiNames{iRoi}]);
end


%%
dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;
i=i+1;
figure(i)
clf
rows=1;
cols=length(roiNames);
for iRoi = 1:length(roiNames)
    subplot(rows,cols,iRoi);
    
    for rwd=1:2
        
        dsErrorsurface(1:trialLength, squeeze(meanResponse(iRoi,rwd,:)), squeeze(stdResponse(iRoi,rwd,:))./sqrt(length(goodSubs)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(squeeze(meanResponse(iRoi,rwd,:)), 'Color', plotColors{rwd}, 'linewidth', 2);
        
    end
    title(roiNames{iRoi})
end

% %% single subject
% i=i+1;
% figure(i)
% clf
% cols=length(subFolders);
% rows=length(roiNames);
% for iRoi = 1:length(roiNames)
%     for iSub = 1:length(subFolders)
%         subplot(rows,cols,(iRoi-1)*cols+iSub);
%         
%         for rwd=1:2
%             plot(squeeze(subResponse(iSub,iRoi,rwd,:)), '.-', 'Color', plotColors{rwd}, 'linewidth', 1,'markersize',20);
%             hold on
%         end
%     end
% end
%%

% if regressGlobalMean
%     i=i+1;
%     figure(i)
%     clf
%     rows=1;
%     cols=length(subFolders);
%     for iSub = 1:length(subFolders)
%         subplot(rows,cols,iSub)
%         for rwd=1:2
%             globalTrial = mean(reshape(globalMean{iSub,rwd},10,[]),2);
%             plot(globalTrial,'Color',plotColors{rwd});
%             hold on
%         end
%         
%     end
% end

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

%%

toc

