close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
toZscore=1;%0 or 1
regressGlobalMean = 0;
ConcatProj=1;
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
dsSurfaceContrast = 0.5;
dsSurfaceAlpha = 0.3;
linewidth = 1;
markersize=10;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
% ROIs = [1 2 7 8];
ROIs = [1 2];
% goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22


%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;


%% bin voxels within each ROI according to eccentricity, then average within bins
%bins by log eccentricity
eccMin = 0.2;
eccMax = 70;
nbins = 12;
% binBorders = (linspace(eccMin^3,eccMax^3,nbins+1))^(1/3)
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
% binBorders = [0 1.3 3 7 20 100];
% binBorders = [0 0.4 1.3 2 3 5 7 12 20 40 100];
nbins = length(binBorders)-1;
% binBorders=exp(linspace(log(eccMin),log(eccMax),nbins+1));
% binBorders = (linspace(eccMin,eccMax,nbins+1));
% minBinVoxels=5;%minimum voxels per bin, otherwise subject is excluded from average
clear subBinTrialResponse subBinResponse binVoxels numBinVoxels meanSubVoxCorr meanVoxCorr
for iSub = 1:length(goodSubs)%length(subdirs)
    for iRoi= ROIs
        for rwd=1:2
            temp = trialCorrectness{goodSubs(iSub),rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);
            temp = trialResponse{goodSubs(iSub),rwd}(:,2:end-1);
            trialResponseVec = temp(:);
            if onlyCorrect ==1 %ONLY CORRECT
                    goodTrials = trialCorrectnessVec==1;
                    %                 numTrials = sum(trialCorrectnessVec);
                elseif onlyCorrect ==2 % ONLY INCORRECT!!!
                    goodTrials = trialCorrectnessVec==0 & trialResponseVec>0;
                else % including all trials with a response
                    goodTrials = trialResponseVec>0;
            end
                
            for ibin=1:nbins
                binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
%                 binVoxels = binVoxels & areas{iSub,iRoi}==1;%ONLY V1
                binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries = nanmean(roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:));%mean timecourse across voxels
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = reshape(binMeanTseries, trialLength, length(binMeanTseries)/trialLength);

%                 if onlyCorrect
%                     subBinTrialResponse{iSub,iRoi,ibin,rwd} = subBinTrialResponse{iSub,iRoi,ibin,rwd}(:,trialCorrectnessVec==1);
%                 end
                
                
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = subBinTrialResponse{iSub,iRoi,ibin,rwd}(:,goodTrials);
                
                
                reshapedTrials = reshape(subBinTrialResponse{iSub,iRoi,ibin,rwd},trialLength,[]);
                if iSub==1
                    allBinTrials{iRoi,ibin,rwd} = reshapedTrials;
                else
                    allBinTrials{iRoi,ibin,rwd} = [allBinTrials{iRoi,ibin,rwd} reshapedTrials];
                end
                %average per subject
                subBinResponse(iSub,iRoi,ibin,rwd,:) = mean(subBinTrialResponse{iSub,iRoi,ibin,rwd},2);
                %trial-by-trial variability per subject
                subBinStd(iSub,iRoi,ibin,rwd,:) = std(subBinTrialResponse{iSub,iRoi,ibin,rwd},0,2);
                
                temp = fft(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                roiBinFftAmp = abs(temp(2,:));
                roiBinFftPh = angle(temp(2,:));
                subBinFftAmpStd(iSub,iRoi,ibin,rwd) = std(roiBinFftAmp);
                subBinFftPhStd(iSub,iRoi,ibin,rwd) = circ_std(roiBinFftPh');
                roiBunFftOther = sum(abs(temp([1 3:end],:)));
                subBinFftOtherStd(iSub,iRoi,ibin,rwd) = std(roiBunFftOther);

            end

        end
    end
end

subBinAmp = squeeze(std(subBinResponse,0,5));%get mean timecourse amplitude per subject
binAmp = squeeze(mean(subBinAmp,1)); %(iRoi,ibin,rwd)
binDiffAmp = squeeze(binAmp(:,:,1) - binAmp(:,:,2));%(iRoi,ibin)

binResponse = squeeze(mean(subBinResponse,1)); % mean timecourse across subjects. binResponse(iRoi,ibin,rwd,timepoint)
binResponseAmp = squeeze(std(binResponse,0,4)); %amplitude of mean timecourse. binResponseAmp(iRoi,ibin,rwd)

%mean trial-by-trial variability
binStdMean = squeeze(mean(subBinStd,1)); % mean variability across subjects. binStd(iRoi,ibin,rwd,timepoint)
binStdStd = squeeze(std(subBinStd,0,1)); % std of variability across subjects. binStd(iRoi,ibin,rwd,timepoint)
binStdDiff = squeeze(binStdMean(:,:,2,:) - binStdMean(:,:,1,:));

%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd pVal_rfx pVal_ffx



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
%         roiTrials{iRoi} = [subTrialResponse{goodSubs(iSub),iRoi,1} subTrialResponse{goodSubs(iSub),iRoi,2}];
        for ibin=1:nbins
            roiBinTrials{iRoi,ibin} = [subBinTrialResponse{iSub,iRoi,ibin,1} subBinTrialResponse{iSub,iRoi,ibin,2}];
        end
    end
    
    for p=1:nperms
        randOrder = randperm(numTrials(iSub,1)+numTrials(iSub,2));
        for iRoi= ROIs%length(roiNames)
            for ibin=1:nbins
                for rwd=1:2
                    permTrials = roiBinTrials{iRoi,ibin}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));% 10 timepoints X #trials
                    permSubMeanTC(iSub,iRoi,ibin,rwd,p,:) = mean(permTrials,2);%average timecourse over trials
                    permSubResp(iSub,iRoi,ibin,rwd,p) = std(squeeze(permSubMeanTC(iSub,iRoi,ibin,rwd,p,:)));%amplitude of task response
                    
                    
                    %trial-by-trial variability per subject
                    permSubStd(iSub,iRoi,ibin,rwd,p,:) = std(permTrials,0,2);%timecourse of variability over trials
                                    
                    %                 permTrials = roiTrials{iRoi}(:,randOrder(firstTrial(rwd):firstTrial(rwd)+numTrials(iSub,rwd)-1));% 10 timepoints X #trials
                    %                 permSubMeanTC(iSub,iRoi,rwd,p,:) = mean(permTrials,2);%average timecourse over trials
                    %                 permSubResp(iSub,iRoi,rwd,p) = std(squeeze(permSubMeanTC(iSub,iRoi,rwd,p,:)));%amplitude of task response
                    
                    %trial-by-trial FFT amp and phase variability per subject
                    temp = fft(permTrials);
                    permFftAmp = abs(temp(2,:));
                    permFftPh = angle(temp(2,:));
                    permSubAmpStd(iSub,iRoi,ibin,rwd,p) = std(permFftAmp);
                    permSubPhStd(iSub,iRoi,ibin,rwd,p) = circ_std(permFftPh');
                    permFftOther = sum(abs(temp([1 3:end],:)));
                    permSubOtherStd(iSub,iRoi,ibin,rwd,p) = std(permFftOther);
                end
            end
        end
    end
    
    for iRoi= ROIs%length(roiNames)
        for ibin=1:nbins
            for rwd = 1:2
                realSubMeanTC(iSub,iRoi,ibin,rwd,:) = mean(subBinTrialResponse{iSub,iRoi,ibin,rwd},2);
                
            end
            realSubDiff(iSub,iRoi,ibin) = std(realSubMeanTC(iSub,iRoi,ibin,1,:)) - std(realSubMeanTC(iSub,iRoi,ibin,2,:));%difference between high std and low std
            permSubDiff(iSub,iRoi,ibin,:) = squeeze(permSubResp(iSub,iRoi,ibin,1,:) - permSubResp(iSub,iRoi,ibin,2,:));%difference between high std and low std
            pVal_sub(iSub, iRoi,ibin) = sum(permSubDiff(iSub,iRoi,ibin,:) > realSubDiff(iSub,iRoi,ibin))/nperms;
            
            permSubStdDiff(iSub,iRoi,ibin,:,:) = squeeze(permSubStd(iSub,iRoi,ibin,2,:,:) - permSubStd(iSub,iRoi,ibin,1,:,:));%difference between high and low reward varability 

            realSubAmpStdDiff(iSub,iRoi,ibin) = subBinFftAmpStd(iSub,iRoi,ibin,2) - subBinFftAmpStd(iSub,iRoi,ibin,1);
            realSubPhStdDiff(iSub,iRoi,ibin) = subBinFftPhStd(iSub,iRoi,ibin,2) - subBinFftPhStd(iSub,iRoi,ibin,1);
            realSubOtherStdDiff(iSub,iRoi,ibin) = subBinFftOtherStd(iSub,iRoi,ibin,2) - subBinFftOtherStd(iSub,iRoi,ibin,1);            

            permSubAmpStdDiff(iSub,iRoi,ibin,:) = squeeze(permSubAmpStd(iSub,iRoi,ibin,2,:) - permSubAmpStd(iSub,iRoi,ibin,1,:));
            permSubPhStdDiff(iSub,iRoi,ibin,:) = squeeze(permSubPhStd(iSub,iRoi,ibin,2,:) - permSubPhStd(iSub,iRoi,ibin,1,:));
            permSubOtherStdDiff(iSub,iRoi,ibin,:) = squeeze(permSubOtherStd(iSub,iRoi,ibin,2,:) - permSubOtherStd(iSub,iRoi,ibin,1,:));

        end

    end
end
permStdDiff = squeeze(mean(permSubStdDiff));%mean over subjects. permStd(iRoi,ibin,p,timepoint)

meanRealAmpStdDiff = squeeze(mean(realSubAmpStdDiff));%mean over subjects
meanRealPhStdDiff = squeeze(mean(realSubPhStdDiff));%mean over subjects
meanRealOtherStdDiff = squeeze(mean(realSubOtherStdDiff));%mean over subjects

meanPermAmpStdDiff = squeeze(mean(permSubAmpStdDiff));%mean over subjects
meanPermPhStdDiff = squeeze(mean(permSubPhStdDiff));%mean over subjects
meanPermOtherStdDiff = squeeze(mean(permSubOtherStdDiff));%mean over subjects

        
hypoth = pVal_sub<0.05;
% pVal_sub(:,ROIs);
for iRoi= ROIs
    for ibin=1:nbins
        %RFX
        meanPermDiff(iRoi,ibin,:) = mean(permSubDiff(:,iRoi,ibin,:));%average over subjects
        meanRealDiff(iRoi,ibin) = mean(realSubDiff(:,iRoi,ibin));%average over subjects
        pVal_rfx(iRoi,ibin) = sum(meanPermDiff(iRoi,ibin,:) > meanRealDiff(iRoi,ibin))/nperms;
        
        %FFX
        for rwd = 1:2
            meanRealTC(iRoi,ibin,rwd,:) = squeeze(mean(realSubMeanTC(:,iRoi,ibin,rwd,:)));%average timecourse over subjects
            meanPermTC(iRoi,ibin,rwd,:,:) = squeeze(mean(permSubMeanTC(:,iRoi,ibin,rwd,:,:)));%average timecourse over subjects
        end
        meanRealDiffStd(iRoi,ibin) = std(squeeze(meanRealTC(iRoi,ibin,1,:))) - std(squeeze(meanRealTC(iRoi,ibin,2,:)));
        for p=1:nperms
            meanPermDiffStd(iRoi,ibin,p) = std(meanPermTC(iRoi,ibin,1,p,:)) - std(meanPermTC(iRoi,ibin,2,p,:));
        end
        pVal_ffx(iRoi,ibin) = sum(meanPermDiffStd(iRoi,ibin,:) > meanRealDiffStd(iRoi,ibin))/nperms;

        %variability
        for t=1:trialLength
            pVal_var_timecourse(iRoi,ibin,t) = sum( permStdDiff(iRoi,ibin,:,t) > binStdDiff(iRoi,ibin,t) )/nperms;
        end
        pVal_var(iRoi,ibin) = sum( mean(permStdDiff(iRoi,ibin,:,:),4) > mean(binStdDiff(iRoi,ibin,:),3) )/nperms;
        
        %FFT variability

        pVal_ampVar(iRoi,ibin) = sum(meanPermAmpStdDiff(iRoi,ibin,:) > meanRealAmpStdDiff(iRoi,ibin))/nperms;
        pVal_phVar(iRoi,ibin) = sum(meanPermPhStdDiff(iRoi,ibin,:) > meanRealPhStdDiff(iRoi,ibin))/nperms;
        pVal_otherVar(iRoi,ibin) = sum(meanPermOtherStdDiff(iRoi,ibin,:) > meanRealOtherStdDiff(iRoi,ibin))/nperms;
    end
end
pVal_rfx
pVal_ffx
pVal_var_timecourse
pVal_var
pVal_ampVar
pVal_phVar
pVal_otherVar
%%
i=1;
figure(i)
clf
rows=2;
cols=length(ROIs);
for iRoi= ROIs
    subplot(rows,cols,iRoi)
    for rwd=1:2
        plot(squeeze(binAmp(iRoi,:,rwd)),'.-','Color', plotColors{rwd},'linewidth',2,'markersize',20);
        hold on
    end
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
    for rwd=1:2
        plot(squeeze(binResponseAmp(iRoi,:,rwd)),'.-', 'Color', plotColors{rwd},'linewidth',2,'markersize',20);
        hold on
    end
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
        for rwd=1:2
            plot(squeeze(binResponse(iRoi,ibin,rwd,:))', 'Color', plotColors{rwd},'linewidth',2);
            hold on
        end
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
i=i+1;
figure(i)
clf
rows=length(ROIs);
cols=nbins;
for iRoi= ROIs
    for ibin=1:nbins
        subplot(rows,cols,ibin + (iRoi-1)*cols)
        for rwd=1:2
            dsErrorsurface(1:trialLength, squeeze(binStdMean(iRoi,ibin,rwd,:)), squeeze(binStdStd(iRoi,ibin,rwd,:))./sqrt(size(subBinStd,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
            hold on
            plot(1:trialLength, squeeze(binStdMean(iRoi,ibin,rwd,:)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
        end
    end
end


% %%
% histbins=1000;
% rows=length(ROIs);
% cols = length(goodSubs);
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
%         
%         
% %%
% i=i+1;
% figure(i)
% clf
% rows=length(ROIs);
% cols = 2;
% for iRoi= ROIs
%     %RFX
%     subplot(rows,cols,1+(iRoi-1)*cols)
%         histogram(squeeze(meanPermDiff(iRoi,:)),histbins); hold all
%         vline(prctile(meanPermDiff(iRoi,:),95),'k');
%         vline(meanRealDiff(iRoi),'r');
%         title(['random effect: ' roiNames{iRoi}]);
%     %FFX
%         subplot(rows,cols,2+(iRoi-1)*cols)
%         histogram(squeeze(meanPermDiffStd(iRoi,:)),histbins); hold all
%         vline(prctile(meanPermDiffStd(iRoi,:),95),'k');
%         vline(meanRealDiffStd(iRoi),'r');
%         title(['fixed effect: ' roiNames{iRoi}]);
% end
% 
% %%
% i=i+1;
% figure(i)
% clf
% plot(squeeze(mean(numBinVoxels,1))','.-','linewidth',2,'markersize',20)
% xticks([0.5:nbins+0.5]);
%     xticklabels(num2str(binBorders','%0.2f'));
%     xlim([0.5 nbins+0.5]);
% title('mean # voxels per eccentricity bin');
% legend(roiNames)
% 
% %%
% clear allEcc
% numAreas=4;
% allEcc = cell(length(ROIs),numAreas);
% for iSub=1:length(goodSubs)
%     for iRoi=ROIs
%        for iarea=0:3
%            allEcc{iRoi,iarea+1} = [allEcc{iRoi,iarea+1}; eccen{iSub,iRoi}(eccen{iSub,iRoi}>0 & areas{iSub,iRoi}==iarea)];
%        end
%     end
% end
% eccHistBins=1000;
% i=i+1;
% figure(i)
% clf
% rows=length(ROIs);
% cols = numAreas;
% for iRoi=ROIs
%     for iarea=0:3
%         subplot(rows,cols,iarea+1+(iRoi-1)*cols);
%         histogram(allEcc{iRoi,iarea+1},eccHistBins);
%         title([roiNames{iRoi} ' area #' num2str(iarea)]);
%         set(gca,'xscale','log')
%     end
% end
%     
% % plot(squeeze(mean(numBinVoxels,1))','.-','linewidth',2,'markersize',20)
% % title('mean # voxels per eccentricity bin');
% % legend(roiNames)
% %%
% dsSurfaceContrast = 0.5;
% dsSurfaceAlpha = 0.3;
% i=i+1;
% figure(i)
% clf
% rows=1;
% cols=length(roiNames);
% for iRoi = 1:length(roiNames)
%     subplot(rows,cols,iRoi);
%     
%     for rwd=1:2
%         
%         dsErrorsurface(1:trialLength, squeeze(meanResponse(rwd,iRoi,:)), squeeze(stdResponse(rwd,iRoi,:))./sqrt(length(goodSubs)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
%         hold on
%         plot(squeeze(meanResponse(rwd,iRoi,:)), 'Color', plotColors{rwd}, 'linewidth', 2);
%         
%     end
%     title(roiNames{iRoi})
% end
% 
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
toc

save([dataFolder 'rwdEcc_concat' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'subFolders', 'roiNames', ...
    'expName','stairThresh','trialLength',...
    'globalMean','regressBetasGlobal','runRwd', ...
    'pVal_rfx', 'pVal_ffx','meanRealDiff','meanRealDiffStd','nperms','binBorders','nbins','numBinVoxels',...
    'pVal_var_timecourse','pVal_var','pVal_ampVar','pVal_phVar','pVal_otherVar');

