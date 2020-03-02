close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

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
TR=1.5;
% trialLength=10;
clear subMeanResponse trialStd meanTrialStd runStd meanRunStd subRwdStd trialFFTamp trialFFTphase meanTrialFFTamp meanTrialFFTphase 
clear stdTrialFFTphase runFFTamp runFFTphase meanRunFFTamp meanRunFFTphase stdRunFFTphase subFFTamp subFFTphase
clear trialMinMax meanTrialMinMax runMinMax meanRunMinMax subRwdMinMax legendCellArray
clear stdRunMax stdTrialMax
clear groupLabels pVal
clear yError yMin yMax permDiff meanPermDiff
clear permSubDiff realSubDiff permSubResp realSubResp meanPermTC permSubMeanTC realSubMeanTC meanRealTC
plotColors = { [1 0 0],[0 0 1],[0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
linewidth = 1;
markersize=10;
fontsize=9;
% ROIs = 1:length(roiNames)-1;
ROIs = 1:length(roiNames);
ROIs = [1 2 7 8];
ROIs = [2];

dsSurfaceContrast = 1;
dsSurfaceAlpha = 0.15;



%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
%we can either get an amplitude for each trial or for each run
ntrials=15;



%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd binCenters


goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
% subjects. 0.0457  0.0140,     0.0553    0.0214

%% bin voxels within each ROI according to eccentricity, then average within bins
%bins by log eccentricity
eccMin = 0.2;
eccMax = 70;
nbins = 12;
% binBorders = (linspace(eccMin^3,eccMax^3,nbins+1))^(1/3)
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);

nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end

clear subBinTrialResponse subBinResponse binVoxels numBinVoxels meanSubVoxCorr meanVoxCorr
for iSub = 1:length(goodSubs)%length(subdirs)
    for iRoi= ROIs
        for rwd=1:2
            temp = trialCorrectness{goodSubs(iSub),rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);
            temp = trialResponse{goodSubs(iSub),rwd}(:,2:end-1);
            trialResponseVec = temp(:);
            temp = trialRT{goodSubs(iSub),rwd}(:,2:end-1);
            trialRTvec = temp(:);
            if onlyCorrect ==1 %ONLY CORRECT
                goodTrials = trialCorrectnessVec==1 & trialRTvec>0 & trialRTvec<maxRT;
                %                 numTrials = sum(trialCorrectnessVec);
            elseif onlyCorrect ==2 % ONLY INCORRECT!!!
                goodTrials = trialCorrectnessVec==0 & trialResponseVec>0 & trialRTvec>0 & trialRTvec<maxRT;
            elseif onlyCorrect == 0 % including all trials with a response
                goodTrials = trialResponseVec>0 & trialRTvec>0 & trialRTvec<maxRT;
            else%onlyCorrect == 4
                goodTrials = ones(size(trialResponseVec));%ALL trials
            end
                
            for ibin=1:nbins
                binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
%                 binVoxels = binVoxels & areas{iSub,iRoi}==1;%ONLY V1
                binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries = nanmean(roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:));%mean timecourse across voxels
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = reshape(binMeanTseries, trialLength, length(binMeanTseries)/trialLength);
                
                
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
                subBinVar(iSub,iRoi,ibin,rwd) = mean(std(subBinTrialResponse{iSub,iRoi,ibin,rwd},0,2));
                
                temp = fft(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                roiBinFftAmp = abs(temp(2,:));
                roiBinFftPh = angle(temp(2,:));
                subBinFftAmpVar(iSub,iRoi,ibin,rwd) = std(roiBinFftAmp);
                subBinFftPhVar(iSub,iRoi,ibin,rwd) = circ_std(roiBinFftPh');
                
                singleTrialStd = std(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                subBinStdVar(iSub,iRoi,ibin,rwd) = std(singleTrialStd);
                subBinStdMean(iSub,iRoi,ibin,rwd) = mean(singleTrialStd);

            end
                %covariation between voxels
                temp = repmat(goodTrials', trialLength,1);
                goodTimepoints = temp(:);
%                 goodTC = roiTC{iSub,iRoi,rwd}.tSeries(binVoxels,goodTimepoints);
                goodTC = roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(:,goodTimepoints);
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

subBinAmp = std(subBinResponse,0,5);%get mean timecourse amplitude per subject.(iRoi,ibin,rwd,1)
binMeanAmp = squeeze(mean(subBinAmp,1)); %mean over subjects. (iRoi,ibin,rwd)
binStdAmp = squeeze(std(subBinAmp,0,1)); %(iRoi,ibin,rwd)
binDiffAmp = squeeze(binMeanAmp(:,:,1) - binMeanAmp(:,:,2));%(iRoi,ibin)

temp = fft(subBinResponse,[],5);
subBinPh = angle(temp(:,:,:,:,2));
binMeanPh = squeeze(circ_mean(subBinPh));
binStdPh = squeeze(circ_std(subBinPh));
subBinPhDiff = squeeze(circ_dist(subBinPh(:,:,:,1),subBinPh(:,:,:,2)));
binMeanPhDiff = squeeze(mean(subBinPhDiff));
binStdPhDiff = squeeze(std(subBinPhDiff));

subBinDiff = squeeze(subBinAmp(:,:,:,1) - subBinAmp(:,:,:,2));%(iSub,iRoi,ibin)
binMeanDiff = squeeze(mean(subBinDiff));
binStdDiff = squeeze(std(subBinDiff));
binResponse = squeeze(mean(subBinResponse,1)); % mean timecourse across subjects. binResponse(iRoi,ibin,rwd,timepoint)
binResponseAmp = squeeze(std(binResponse,0,4)); %amplitude of mean timecourse. binResponseAmp(iRoi,ibin,rwd)



binVarMean = squeeze(mean(subBinVar,1)); % mean variability across subjects. binStdMean(iRoi,ibin,rwd,timepoint)
binVarStd = squeeze(std(subBinVar,0,1)); % std of variability across subjects. binStdStd(iRoi,ibin,rwd,timepoint)
binSubVarDiff = squeeze(subBinVar(:,:,:,2) - subBinVar(:,:,:,1));
binVarDiffMean = squeeze(mean(binSubVarDiff));
binVarDiffStd = squeeze(std(binSubVarDiff));

binStdVarMean = squeeze(mean(subBinStdVar));%mean trial std variability
binStdVarStd = squeeze(std(subBinStdVar));%std across subjects, of trial std variability
binSubStdVarDiff = squeeze(subBinStdVar(:,:,:,2) - subBinStdVar(:,:,:,1));
binStdVarDiffMean = squeeze(mean(binSubStdVarDiff));
binStdVarDiffStd = squeeze(std(binSubStdVarDiff));

% binStdMean = squeeze(mean(subBinStdMean));%mean trial std mean

binFftAmpVarMean = squeeze(mean(subBinFftAmpVar,1)); % mean variability across subjects. binStdMean(iRoi,ibin,rwd,timepoint)
binFftAmpVarStd = squeeze(std(subBinFftAmpVar,0,1)); % std of variability across subjects. binStdStd(iRoi,ibin,rwd,timepoint)
binSubFftAmpVarDiff = squeeze(subBinFftAmpVar(:,:,:,2) - subBinFftAmpVar(:,:,:,1));
binFftAmpVarDiffMean = squeeze(mean(binSubFftAmpVarDiff));
binFftAmpVarDiffStd = squeeze(std(binSubFftAmpVarDiff));

binFftPhVarMean = squeeze(mean(subBinFftPhVar,1)); % mean variability across subjects. binStdMean(iRoi,ibin,rwd,timepoint)
binFftPhVarStd = squeeze(std(subBinFftPhVar,0,1)); % std of variability across subjects. binStdStd(iRoi,ibin,rwd,timepoint)
binSubFftPhVarDiff = squeeze(subBinFftPhVar(:,:,:,2) - subBinFftPhVar(:,:,:,1));
binFftPhVarDiffMean = squeeze(mean(binSubFftPhVarDiff));
binFftPhVarDiffStd = squeeze(std(binSubFftPhVarDiff));




%%
i=0;

rows=2;
cols = 9;
subplots = {1:3, 4:6, 7 , 9, cols+1:cols+3, cols+4:cols+6, cols+7 , cols+9};

%% bar plot of subjects' amplitudes
lineLength = 0.2;
lineWidth = 2;
markerSize = 10;





%%
rows=3;
cols = 9;
subplots = {1:3, 4:6, 7 , 9};
% iRoi = ROIs(1);

for r = 1:length(ROIs)
    iRoi = ROIs(r);
    i=i+1;
    figure(i)
    clf
    
    % timepoint variability
    subplot(rows,cols,subplots{1});
%     for rwd=1:2
%         dsErrorsurface(binCenters, squeeze(binVarMean(iRoi,:,rwd)), squeeze(binVarStd(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
%         hold on
%     end
%     for rwd=1:2
%         plot(binCenters, squeeze(binVarMean(iRoi,:,rwd)),'.','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
%     end
%     ylabel('mean variability (std)');
    
    %%
    %set canonical HRF
    upsampleFactor = 10;
    TR = 1.5;
    modelParams = struct;
    sampleDuration = TR/upsampleFactor;
    sampleDelay=sampleDuration/2;
    defaultParams=1;
%     modelParams.x = 6;%9;
%     modelParams.y = 16;
%     modelParams.z = 6;
%      modelParams.x = 12;
%     modelParams.y = 16;
%     modelParams.z = 6;

    modelParams.x = 6;%9;
modelParams.y = 12;
modelParams.z = 100;

    [modelParams hrfModel] = hrfDoubleGamma(modelParams,sampleDuration,sampleDelay,defaultParams);
    T=200;
    initT=1;
    simResp = hrfModel(initT:initT+T-1);
    linewidth = 1;
    linecolor = 0.6*[1 1 1];
%     tickcolor = 0.2*[1 1 1];
    tickcolor = 0.6*[0.0 1 0.9];
    nCurves = 5;
    minAmp=0.4;
    amplitudes = linspace(minAmp, 0.9,nCurves);
    minDelay = 0;
    maxDelay = 0.1;
    delays = linspace(minDelay,maxDelay,nCurves);
%     delays = delays(randperm(length(delays)));
    delays = delays([2 4 1 5 3]);
%     amplitudes = [1 0.8 0.6];
%     delays = [0 0.3 0.6];
    clear curves
    for icurve = 1:nCurves
       curves(icurve,:) = amplitudes(icurve) *  circshift(simResp, ceil(delays(icurve)*T));
    end
    
    [maxAmp maxT] = max(curves');
    Tscale = 1/50;
    ampScale = 1/25;
%     horizX = [(T*Tscale)*ones(size(maxT)); (T*Tscale)*ones(size(maxT)) + T*Tscale];
    horizX = [zeros(size(maxT)); T*Tscale*ones(size(maxT))];
    horizXlong = [zeros(size(maxT));  maxT];
    horizY = repmat(maxAmp,[2 1]);


    vertX = repmat(maxT, [2 1]);
    temp1 = (-max(maxAmp)*ampScale)*ones(size(maxT));
    temp2 = (-max(maxAmp)*ampScale)*ones(size(maxT)) -max(maxAmp)*ampScale;
    vertY = [zeros(size(maxT)); max(maxAmp)*ampScale*ones(size(maxT))];
    vertYlong = [zeros(size(maxT)); maxAmp];
    
    N = 15;

    vertXrand = (maxDelay-minDelay)*0.3*T*repmat(randn(1,N), [2 1]) + mean(maxT);
    temp1 = (-max(maxAmp)*ampScale)*ones(1,N);
    temp2 = (-max(maxAmp)*ampScale)*ones(1,N) -max(maxAmp)*ampScale;
    vertYrand = [zeros(1,N); max(maxAmp)*ampScale*ones(1,N)];
    
    horizXrand = [zeros(1,N); (T*Tscale)*ones(1,N)];
    horizYrand = minAmp*0.3*max(maxAmp)*repmat(randn(1,N), [2 1]) + mean(maxAmp);
    xmin = 0;
    xmax = (maxDelay+0.4)*T;
    xmax = 1+ceil(1+ xmax/(trialLength-1))*(trialLength-1);
    for irow=1:3
        subplot(rows,cols,(irow-1)*cols+subplots{1});
        cla
        plot(curves','Color',linecolor,'linewidth',linewidth);
        hold on
        %     line(horizX,horizY,'color',0.2*[1 1 1],'linewidth',linewidth);
        %     line(horizXlong,horizY,'linestyle','--','color',0.2*[1 1 1],'linewidth',linewidth/3);
        
%         xlim([(minDelay+0.1)*T (maxDelay+0.3)*T]);

    end
    
    subplot(rows,cols,subplots{1});
    meanCurves = mean(curves);
    stdCurves = std(curves);
    stdXvals = linspace(xmin+1,xmax,trialLength);
    errorbar(stdXvals, meanCurves(stdXvals), stdCurves(stdXvals),'color',tickcolor,'linewidth',linewidth,'linestyle',':'); 
    
    subplot(rows,cols,cols+subplots{1});
    line(horizXrand,horizYrand,'linestyle','-','color',tickcolor,'linewidth',linewidth*2);
    line(horizXlong,horizY,'linestyle',':','color',tickcolor,'linewidth',linewidth);
    line(horizX,horizY,'color',tickcolor,'linewidth',linewidth*2);
    subplot(rows,cols,2*cols+subplots{1});
    line(vertXrand,vertYrand,'linestyle','-','color',tickcolor,'linewidth',linewidth*2);
    line(vertX,vertYlong,'linestyle',':','color',tickcolor,'linewidth',linewidth);
    line(vertX,vertY,'color',tickcolor,'linewidth',linewidth*2);


    
%     ylim([0 1.05*max(maxAmp)]);
    %%
    
    
    subplot(rows,cols,subplots{2});
    dsErrorsurface(binCenters, binVarDiffMean(iRoi,:), binVarDiffStd(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
    hold on
    plot(binCenters, binVarDiffMean(iRoi,:),'k.','linewidth',linewidth,'markersize',markersize);
%     ylabel('\Delta timepoint variability (std)');
ylabel('high reward - low reward');
hline(0);
    
%     % STD amplitude variability
%     subplot(rows,cols,cols+subplots{1});
%     for rwd=1:2
%         dsErrorsurface(binCenters, squeeze(binStdVarMean(iRoi,:,rwd)), squeeze(binStdVarStd(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
%         hold on
%     end
%     for rwd=1:2
%         plot(binCenters, squeeze(binStdVarMean(iRoi,:,rwd)),'.','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
%     end
%     ylabel('amplitude variability (std)');
    subplot(rows,cols,cols+subplots{2});
    dsErrorsurface(binCenters, binStdVarDiffMean(iRoi,:), binStdVarDiffStd(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
    hold on
    plot(binCenters, binStdVarDiffMean(iRoi,:),'k.','linewidth',linewidth,'markersize',markersize);
%     ylabel('\Delta amplitude variability (std)');
ylabel('high reward - low reward');
    hline(0);
    % FFT amplitude variability
    % subplot(rows,cols,cols+subplots{1});
    % for rwd=1:2
    %     dsErrorsurface(binCenters, squeeze(binFftAmpVarMean(iRoi,:,rwd)), squeeze(binFftAmpVarStd(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    %     hold on
    %     plot(binCenters, squeeze(binFftAmpVarMean(iRoi,:,rwd)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
    % end
    % ylabel('amplitude variability (std)');
    % subplot(rows,cols,cols+subplots{2});
    % dsErrorsurface(binCenters, binFftAmpVarDiffMean(iRoi,:), binFftAmpVarDiffStd(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
    % hold on
    % plot(binCenters, binFftAmpVarDiffMean(iRoi,:),'k.-','linewidth',linewidth,'markersize',markersize);
    % ylabel('\Delta amplitude variability (std)');
    % hline(0);
    
    % FFT phase variability
    
%     subplot(rows,cols,2*cols+subplots{1});
%     for rwd=1:2
%         dsErrorsurface(binCenters, squeeze(binFftPhVarMean(iRoi,:,rwd)), squeeze(binFftPhVarStd(iRoi,:,rwd))./sqrt(size(subBinAmp,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
%         hold on
%     end
%     for rwd=1:2
%         plot(binCenters, squeeze(binFftPhVarMean(iRoi,:,rwd)),'.','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
%     end
%     ylabel('temporal variability (std)');
    subplot(rows,cols,2*cols+subplots{2});
    dsErrorsurface(binCenters, binFftPhVarDiffMean(iRoi,:), binFftPhVarDiffStd(iRoi,:)./sqrt(size(subBinAmp,1)), [0 0 0],dsSurfaceAlpha);
    hold on
    plot(binCenters, binFftPhVarDiffMean(iRoi,:),'k.','linewidth',linewidth,'markersize',markersize);
%     ylabel('\Delta temporal variability (std)');
ylabel('high reward - low reward');
    hline(0);
    
    
    
    % get subject colors for bar plot
    clear smallSubAmp
    for iSub=1:length(goodSubs)
        for rwd=1:2
            smallSubAmp(iSub,rwd) = std(squeeze(subResponse(goodSubs(iSub),iRoi,rwd,:)));
        end
    end
    smallSubDiff = smallSubAmp(:,1) - smallSubAmp(:,2);
    minSubDiff = min(smallSubDiff);
    maxSubDiff= max(smallSubDiff);
    subjects = size(smallSubAmp,1);
    rewards = size(smallSubAmp,2);
    [rwdNum, subNum] = meshgrid(1:rewards, 1:subjects);
    scatterCmap=cool;
    %     scatterCmap = scatterCmap(end:-1:1,:);%invert color map
    colormap(scatterCmap);
    l = size(scatterCmap,1);
    for iSub=1:length(goodSubs)
        subColor(iSub,:) = scatterCmap(1+floor((smallSubDiff(iSub) - minSubDiff)*(l-1)/(maxSubDiff-minSubDiff)),:);
    end
    
    % compute subject variability for whole ROI
    clear subVar subMeanVar roiBinFftAmp roiBinFftPh roiBinFftOther subBinFftAmpVar subBinFftPhVar subBinFftOtherVar
    for iSub = 1:length(goodSubs)%length(subdirs)
        for rwd=1:2
            reshapedTrials = reshape(subTrialResponse{goodSubs(iSub),iRoi,rwd},10,[]);
            %trial-by-trial variability per subject
            subVar(iSub,rwd,:) = std(reshapedTrials,0,2);
            subMeanVar(iSub,rwd) = mean(squeeze(subVar(iSub,rwd,:)));
            temp = fft(reshapedTrials);
            roiBinFftAmp = abs(temp(2,:));
            roiBinFftPh = angle(temp(2,:));
            roiBinFftOther = sum(abs(temp([1 3:end],:)));
            subBinFftAmpVar(iSub,rwd) = std(roiBinFftAmp);
            subBinFftPhVar(iSub,rwd) = circ_std(roiBinFftPh');
            subBinFftOtherVar(iSub,rwd) = std(roiBinFftOther);
            
            singleTrialStd = std(reshapedTrials);
            subStdVar(iSub,rwd) = std(singleTrialStd);
            subStdMean(iSub,rwd) = mean(singleTrialStd);
        end
    end
    
    % bar plot of subjects' timepoint variability
    lineLength = 0.2;
    lineWidth = 2;
    markerSize = 10;
    
    subplot(rows,cols,subplots{3})
    for iSub=1:length(goodSubs)
        plot(1:rewards,subMeanVar(iSub,:),'Color',subColor(iSub,:),'linewidth',linewidth);
        hold on
    end
    scatter(rwdNum(:),subMeanVar(:),markerSize, [subColor; subColor] ,'filled');
    % MEAN
    for iRwd=1:rewards
        scatter(iRwd,mean(subMeanVar(:,iRwd)),markerSize*4,[0 0 0],'filled');
        plot(1:rewards,mean(subMeanVar),'Color',[0 0 0],'linewidth', 2*linewidth);
    end
    ylabel('timepoint variability (std)');
    
    % bar plot of subjects' STD Amplitude variability
    subplot(rows,cols,cols+subplots{3})
    for iSub=1:length(goodSubs)
        plot(1:rewards,subStdVar(iSub,:),'Color',subColor(iSub,:),'linewidth',linewidth);
        hold on
    end
    scatter(rwdNum(:),subStdVar(:),markerSize, [subColor; subColor] ,'filled');
    % MEAN
    for iRwd=1:rewards
        scatter(iRwd,mean(subStdVar(:,iRwd)),markerSize*4,[0 0 0],'filled');
        plot(1:rewards,mean(subStdVar),'Color',[0 0 0],'linewidth', 2*linewidth);
    end
    ylabel('amplitude variability (std)');
    
    % bar plot of subjects' FFT Amplitude variability
    % subplot(rows,cols,cols+subplots{3})
    % for i=1:length(goodSubs)
    %     plot(1:rewards,subBinFftAmpVar(i,:),'Color',subColor(i,:),'linewidth',linewidth);
    %     hold on
    % end
    % scatter(rwdNum(:),subBinFftAmpVar(:),markerSize, [subColor; subColor] ,'filled');
    % % MEAN
    % for i=1:rewards
    %     scatter(i,mean(subBinFftAmpVar(:,i)),markerSize*4,[0 0 0],'filled','d');
    %     plot(1:rewards,mean(subBinFftAmpVar),'Color',[0 0 0],'linewidth', 2*linewidth);
    % end
    % ylabel('amplitude variability (std)');
    
    % bar plot of subjects' FFT Phase variability
    subplot(rows,cols,2*cols+subplots{3})
    for iSub=1:length(goodSubs)
        plot(1:rewards,subBinFftPhVar(iSub,:),'Color',subColor(iSub,:),'linewidth',linewidth);
        hold on
    end
    scatter(rwdNum(:),subBinFftPhVar(:),markerSize, [subColor; subColor] ,'filled');
    % MEAN
    for iRwd=1:rewards
        scatter(iRwd,mean(subBinFftPhVar(:,iRwd)),markerSize*4,[0 0 0],'filled');
        plot(1:rewards,mean(subBinFftPhVar),'Color',[0 0 0],'linewidth', 2*linewidth);
    end
    ylabel('temporal variability (std)');
    
    
    
  %%  
    % formatting
    for r=1:rows
        for isubplot=1:length(subplots)-1
            subplot(rows,cols,(r-1)*cols+subplots{isubplot})
            %         ylabel('mean variability (std)');

            if isubplot<3
                if isubplot==1
                    
                    
%                     xlim([xmin xmax]);
%                     ylim([0 max(maxAmp)]);
                    yticklabel = {'0','1'};
                    ytick = [0 max(maxAmp)*1.2];
                    xticklabel = {'0','15'};
                    xtick = [0 xmax];
                    drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -8/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
                    'labelFontSize',fontsize,'yTickLabel',yticklabel,'yTick',ytick,'xTickLabel',xticklabel,'xTick',xtick,...
                    'yLabel','fMRI response','xLabel','time (sec)');
%                     ylabel('fMRI response');
%                     xlabel('time (sec)');
                else
%                     xlabel('eccentricity (deg)');
                     drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
                    'labelFontSize',fontsize,'xLabel','eccentricity (deg)');
                    
                end
                
                axis square
            else
%                 xlabel('reward');
                drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -70/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
                    'xTick',[1 2], 'xTickLabel', {'high','low'},...
                    'xAxisMin', 1 ,'xAxisMax', 2,'yAxisOffset',-0.5,'labelFontSize',fontsize,...
                    'yAxisMajorTickLen',-4/32,'xLabel','reward');
            end
        end
    end
    
    
    set(gcf,'position',[10 10 26 	22]);
    
    print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig5_' ConcatProjStr '.pdf']);
end

colorSubFolders = subFolders;
colorGoodSubs = goodSubs;
save(fullfile(dataFolder,'subColors.mat'),'subColor', 'scatterCmap', 'colorSubFolders', 'colorGoodSubs');
%%

% [rmat pmat] = corr([smallSubDiff,subMeanVar(:,2)-subMeanVar(:,1),subBinFftAmpVar(:,2)-subBinFftAmpVar(:,1),...
%     subBinFftPhVar(:,2)-subBinFftPhVar(:,1), subStdMean(:,2)-subStdMean(:,1),subStdVar(:,2)-subStdVar(:,1)]);
% lbls = {'std amp','timepoint var','FFT amp var','FFT ph var','std amp mean','std amp var'};
% 
% pmat


