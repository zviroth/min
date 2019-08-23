onlyCorrect=2;%1=correct,2=incorrect,0=all
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
load([dataFolder 'rwdTC_concat' onlyCorrectString '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', 'meanResponse', 'stdResponse',...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength');
allTrials
keyboard
% trialLength=10;
clear subMeanResponse trialStd meanTrialStd runStd meanRunStd subRwdStd trialFFTamp trialFFTphase meanTrialFFTamp meanTrialFFTphase 
clear stdTrialFFTphase runFFTamp runFFTphase meanRunFFTamp meanRunFFTphase stdRunFFTphase subFFTamp subFFTphase
clear trialMinMax meanTrialMinMax runMinMax meanRunMinMax subRwdMinMax legendCellArray
clear stdRunMax stdTrialMax
clear groupLabels pVal
clear yError yMin yMax permDiff meanPermDiff
clear permSubDiff realSubDiff permSubResp realSubResp meanPermTC permSubMeanTC realSubMeanTC meanRealTC
plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};

% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
ROIs = [1 2 7 8];
ROIs = [1 2];



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
goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% subjects. 0.0457  0.0140,     0.0553    0.0214

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




%% PERMUTATIONS
tic

for iSub = 1:length(goodSubs)%length(subdirs)
    for rwd=1:2

        numTrials(iSub,rwd) = size(subTrialResponse{goodSubs(iSub),1,rwd},2);%may be different number of trials for low and high reward!
        if rwd==1
            firstTrial(rwd)=1;
        else %rwd==2
            firstTrial(rwd)=numTrials(iSub,1)+1;
        end
    end
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
        

    end
end
hypoth = pVal_sub<0.05;
% pVal_sub(:,ROIs);
for iRoi= ROIs
    %RFX
    meanPermDiff(iRoi,:) = mean(permSubDiff(:,iRoi,:));%average over subjects
    meanRealDiff(iRoi) = mean(realSubDiff(:,iRoi));%average over subjects
    pVal_rfx(iRoi) = sum(meanPermDiff(iRoi,:) > meanRealDiff(iRoi))/nperms;
    
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
end
pVal_rfx
pVal_ffx

%% bin voxels within each ROI according to eccentricity, then average within bins
%bins by log eccentricity
eccMin = 0.2;
eccMax = 70;
nbins = 5;
% binBorders = (linspace(eccMin^3,eccMax^3,nbins+1))^(1/3)
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
binBorders = [0 1.3 3 7 20 100];
binBorders = [0 0.4 1.3 2 3 5 7 12 20 40 100];
nbins = length(binBorders)-1;
% binBorders=exp(linspace(log(eccMin),log(eccMax),nbins+1));
% binBorders = (linspace(eccMin,eccMax,nbins+1));
% minBinVoxels=5;%minimum voxels per bin, otherwise subject is excluded from average
clear subBinTrialResponse subBinResponse binVoxels numBinVoxels meanSubVoxCorr meanVoxCorr
for iSub = 1:length(goodSubs)%length(subdirs)
    for iRoi= ROIs
        for rwd=1:2
            temp = trialCorrectness{iSub,rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);
            temp = trialResponse{iSub,rwd}(:,2:end-1);
            trialResponseVec = temp(:);
            
            for ibin=1:nbins
                binVoxels = eccen{iSub,iRoi}>binBorders(ibin) & eccen{iSub,iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
%                 binVoxels = binVoxels & areas{iSub,iRoi}==1;%ONLY V1
                binVoxels = binVoxels & areas{iSub,iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries = nanmean(roiTC{iSub,iRoi,rwd}.tSeries(binVoxels,:));%mean timecourse across voxels
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = reshape(binMeanTseries, trialLength, length(binMeanTseries)/trialLength);

%                 if onlyCorrect
%                     subBinTrialResponse{iSub,iRoi,ibin,rwd} = subBinTrialResponse{iSub,iRoi,ibin,rwd}(:,trialCorrectnessVec==1);
%                 end
                
                if onlyCorrect ==1 %ONLY CORRECT
                    goodTrials = trialCorrectnessVec==1;
                    %                 numTrials = sum(trialCorrectnessVec);
                elseif onlyCorrect ==2 % ONLY INCORRECT!!!
                    goodTrials = trialCorrectnessVec==0 & trialResponseVec>0;
                else % including all trials with a response
                    goodTrials = trialResponseVec>0;
                end
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = subBinTrialResponse{iSub,iRoi,ibin,rwd}(:,goodTrials);
                
                
                reshapedTrials = reshape(subBinTrialResponse{iSub,iRoi,ibin,rwd},trialLength,[]);
                if iSub==1
                    allBinTrials{iRoi,ibin,rwd} = reshapedTrials;
                else
                    allBinTrials{iRoi,ibin,rwd} = [allBinTrials{iRoi,ibin,rwd} reshapedTrials];
                end
                %average per subject
                subBinResponse(iSub,iRoi,ibin,rwd,:) = mean(subBinTrialResponse{iSub,iRoi,ibin,rwd},2);
                
                

            end
                %covariation between voxels
                temp = repmat(goodTrials', trialLength,1);
                goodTimepoints = temp(:);
%                 goodTC = roiTC{iSub,iRoi,rwd}.tSeries(binVoxels,goodTimepoints);
                goodTC = roiTC{iSub,iRoi,rwd}.tSeries(:,goodTimepoints);
                %   maybe remove mean trial per voxel?
                temp=triu(corr(goodTC'));
                meanSubVoxCorr(iSub,iRoi,rwd) = nanmean(temp(:));
        end
    end
end
meanVoxCorr = squeeze(nanmean(meanSubVoxCorr));
fisherZ = atan(meanSubVoxCorr);
for iRoi= ROIs
    [h,  roiCorrPval(iRoi)] = ttest(fisherZ(:,iRoi,1)-fisherZ(:,iRoi,2));

end

subBinAmp = squeeze(std(subBinResponse,0,5));%get mean timecourse amplitude per subject
binAmp = squeeze(mean(subBinAmp,1)); %(iRoi,ibin,rwd)
binDiffAmp = squeeze(binAmp(:,:,1) - binAmp(:,:,2));%(iRoi,ibin)

binResponse = squeeze(mean(subBinResponse,1)); % mean timecourse across subjects. binResponse(iRoi,ibin,rwd,timepoint)
binResponseAmp = squeeze(std(binResponse,0,4)); %amplitude of mean timecourse. binResponseAmp(iRoi,ibin,rwd)


%%
i=1;
figure(i)
clf
rows=2;
cols=length(ROIs);
for iRoi= ROIs
    subplot(rows,cols,iRoi)
    plot(squeeze(binAmp(iRoi,:,:)),'.-','linewidth',2,'markersize',20);
    xticks([0.5:nbins+0.5]);
    xticklabels(num2str(binBorders','%0.2f'));
    xlim([0.5 nbins+0.5]);
    title(['mean amp: ' roiNames{iRoi}])
%     set(gca,'xscale','log')

    subplot(rows,cols,iRoi+cols)
    plot(squeeze(binAmp(iRoi,:,1) -binAmp(iRoi,:,2)),'k.-','linewidth',2,'markersize',20);
    xticks([0.5:nbins+0.5]);
    xticklabels(num2str(binBorders','%0.2f'));
    xlim([0.5 nbins+0.5]);
    title(['H minus L'])
end
% legend('H','L')

%%
i=i+1;
figure(i)
clf
rows=2;
cols=length(ROIs);
for iRoi= ROIs
    subplot(rows,cols,iRoi)
    plot(squeeze(binResponseAmp(iRoi,:,:)),'.-','linewidth',2,'markersize',20);
    xticks([0.5:nbins+0.5]);
    xticklabels(num2str(binBorders','%0.2f'));
    xlim([0.5 nbins+0.5]);
    title(['amp of mean: ' roiNames{iRoi}])
%     set(gca,'xscale','log')

    subplot(rows,cols,iRoi+cols)
    plot(squeeze(binResponseAmp(iRoi,:,1) -binResponseAmp(iRoi,:,2)),'k.-','linewidth',2,'markersize',20);
    xticks([0.5:nbins+0.5]);
    xticklabels(num2str(binBorders','%0.2f'));
    xlim([0.5 nbins+0.5]);
    title(['H minus L'])
end
legend('H','L')


%%
i=i+1;
figure(i)
clf
rows=length(ROIs);
cols=nbins;
for iRoi= ROIs
    for ibin=1:nbins
        subplot(rows,cols,ibin + (iRoi-1)*cols)
        plot(squeeze(binResponse(iRoi,ibin,:,:))','linewidth',2);
        title(['mean timecourse: ' roiNames{iRoi} ' ecc ' num2str(ibin)])
    end
%     subplot(rows,cols,iRoi)
%     plot(squeeze(binResponseAmp(iRoi,:,:)),'.-','linewidth',2,'markersize',20);
%     xticks([0.5:nbins+0.5]);
%     xticklabels(num2str(binBorders','%0.2f'));
%     xlim([0.5 nbins+0.5]);
%     title(['amp of mean: ' roiNames{iRoi}])
%     set(gca,'xscale','log')
end
legend('H','L')




%%
histbins=1000;
rows=length(ROIs);
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
rows=length(ROIs);
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
i=i+1;
figure(i)
clf
plot(squeeze(mean(numBinVoxels,1))','.-','linewidth',2,'markersize',20)
xticks([0.5:nbins+0.5]);
    xticklabels(num2str(binBorders','%0.2f'));
    xlim([0.5 nbins+0.5]);
title('mean # voxels per eccentricity bin');
legend(roiNames)

%%
clear allEcc
numAreas=4;
allEcc = cell(length(ROIs),numAreas);
for iSub=1:length(goodSubs)
    for iRoi=ROIs
       for iarea=0:3
           allEcc{iRoi,iarea+1} = [allEcc{iRoi,iarea+1}; eccen{iSub,iRoi}(eccen{iSub,iRoi}>0 & areas{iSub,iRoi}==iarea)];
       end
    end
end
eccHistBins=1000;
i=i+1;
figure(i)
clf
rows=length(ROIs);
cols = numAreas;
for iRoi=ROIs
    for iarea=0:3
        subplot(rows,cols,iarea+1+(iRoi-1)*cols);
        histogram(allEcc{iRoi,iarea+1},eccHistBins);
        title([roiNames{iRoi} ' area #' num2str(iarea)]);
        set(gca,'xscale','log')
    end
end
    
% plot(squeeze(mean(numBinVoxels,1))','.-','linewidth',2,'markersize',20)
% title('mean # voxels per eccentricity bin');
% legend(roiNames)
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
        
        dsErrorsurface(1:trialLength, squeeze(meanResponse(rwd,iRoi,:)), squeeze(stdResponse(rwd,iRoi,:))./sqrt(length(goodSubs)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(squeeze(meanResponse(rwd,iRoi,:)), 'Color', plotColors{rwd}, 'linewidth', 2);
        
    end
    title(roiNames{iRoi})
end


%%
toc



