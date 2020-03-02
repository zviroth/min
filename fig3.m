close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean=0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

randSeed = rng;
load(fullfile(dataFolder,'randSeed.mat'),'randSeed');
TR=1.5;
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
fmriUnits = '% change image intensity';
if toZscore
    zScoreString = '_zscored';
    fmriUnits = 'std image intensity';
end
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
dsSurfaceContrast = 1;
dsSurfaceAlpha = 0.15;
linewidth = 1;
markersize=10;
fontsize = 9;
ROIs = 1:length(roiNames);
ROIs = [1 2];
% goodSubs = 1:length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22

ntrials=15;


%% bin voxels within each ROI according to eccentricity, then average within bins
%bins by log eccentricity
% eccMin = 0.2;
% eccMax = 70;
% nbins = 12;
% binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
binBorders = [0 1 10 100];
nbins = length(binBorders)-1;
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
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = subBinTrialResponse{iSub,iRoi,ibin,rwd}(:,goodTrials);
                
                subBinRuns = reshape(binMeanTseries, trialLength*ntrials, []);
                subBinRunMean(iSub,iRoi,ibin,rwd,:) = mean(subBinRuns,2);
                subBinRunStd(iSub,iRoi,ibin,rwd,:) = std(subBinRuns,0,2);
                subBinFft(iSub,iRoi,ibin,rwd,:) = fft(squeeze(subBinRunMean(iSub,iRoi,ibin,rwd,:)));
                
                
                % REMOVE AMPLITUDE VARIABILITY
                reshapedTrials = reshape(subBinTrialResponse{iSub,iRoi,ibin,rwd},trialLength,[]);
                temp = fft(reshapedTrials);
                roiBinFftAmp = abs(temp(2,:));
                roiBinFftBaseline = abs(temp(1,:));
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = reshapedTrials;
                
                
                
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
                
                singleTrialStd = std(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                subBinStdVar(iSub,iRoi,ibin,rwd) = std(singleTrialStd);
                subBinStdMean(iSub,iRoi,ibin,rwd) = mean(singleTrialStd);

            end

        end
    end
end

subBinAmp = squeeze(std(subBinResponse,0,5));%get mean timecourse amplitude per subject
binAmp = squeeze(mean(subBinAmp,1)); %(iRoi,ibin,rwd)
binDiffAmp = squeeze(binAmp(:,:,1) - binAmp(:,:,2));%(iRoi,ibin)

binResponse = squeeze(mean(subBinResponse,1)); % mean timecourse across subjects. binResponse(iRoi,ibin,rwd,timepoint)
binResponseStd = squeeze(std(subBinResponse,1)); % mean timecourse across subjects. binResponse(iRoi,ibin,rwd,timepoint)
binResponseAmp = squeeze(std(binResponse,0,4)); %amplitude of mean timecourse. binResponseAmp(iRoi,ibin,rwd)

%mean trial-by-trial variability
binStdMean = squeeze(mean(subBinStd,1)); % mean variability across subjects. binStd(iRoi,ibin,rwd,timepoint)
binStdStd = squeeze(std(subBinStd,0,1)); % std of variability across subjects. binStd(iRoi,ibin,rwd,timepoint)
binStdDiff = squeeze(binStdMean(:,:,2,:) - binStdMean(:,:,1,:));

%% PERMUTATIONS
% first combine all trials across all subjects
% we are recreating allTrials, to include only good subjects
clear allTrials meanRwdStd pVal_rfx pVal_ffx

tic


  

%% PLOT timecourse and FFT for s0008
i=i+1; figure(i); clf
iSub=2;
cols=nbins;
rows=2;
iRoi=2;
%timecourse
for ibin=1:nbins
    subplot(rows,cols,ibin)
    for rwd=1:2
        dsErrorsurface(TR*(0:trialLength-1), squeeze(subBinResponse(iSub,iRoi,ibin,rwd,:)), squeeze(subBinStd(iSub,iRoi,ibin,rwd,:))./sqrt(size(subBinTrialResponse{iSub,iRoi,ibin,rwd},2)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(TR*(0:trialLength-1), squeeze(subBinResponse(iSub,iRoi,ibin,rwd,:)),'-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
    end
    xlabel('time (sec)');
        t = {'fMRI response'; ['(' fmriUnits ')'] };
        if ibin==1
       ylabel(t);
        end
%     ylabel(['fMRI response (' fmriUnits ')']);
end
%FFT spectrum
TR=1.5;
nframes = size(subBinFft,5);
runFFT = abs(subBinFft);
runFFT = runFFT / (nframes/2);
frequencies = [0:nframes-1]/(nframes*TR);
runFFT = runFFT(:,:,:,:,1:floor(nframes/2));
frequencies = frequencies(1:floor(nframes/2));

for ibin=1:nbins
    subplot(rows,cols,ibin+cols)
%     subBinRunMean(iSub,iRoi,ibin,rwd,:)
%     subBinRunFft(iSub,iRoi,ibin,rwd,:)
    for rwd=1:2
        plot(1:(trialLength*ntrials)/2, squeeze(squeeze(runFFT(iSub,iRoi,ibin,rwd,:))),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
        hold on
    end
    xlabel('frequency (cycles/run)');
            t = {'Fourier amplitude'; ['(' fmriUnits ')'] };
            if ibin==1
                ylabel(t);
            end
            %     ylabel(['Fourier amplitude  (' fmriUnits ')']);
end

%% Plot timecourse and FFT for group mean
binFft = squeeze(mean(runFFT));
% binFftStd = squeeze(std(runFFT));


i=i+1; figure(i); clf
for ibin=1:nbins
    subplot(rows,cols,ibin)
    for rwd=1:2
        dsErrorsurface(TR*(0:trialLength-1), squeeze(binResponse(iRoi,ibin,rwd,:)), squeeze(binResponseStd(iRoi,ibin,rwd,:))./sqrt(size(subBinResponse,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
        hold on
        plot(TR*(0:trialLength-1), squeeze(binResponse(iRoi,ibin,rwd,:)),'-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
    end
    xlabel('time (sec)');
    t = {'fMRI response'; ['(' fmriUnits ')'] };
    if ibin==1
       ylabel(t);
    end
%     ylabel(['fMRI response (' fmriUnits ')']);
end

for ibin=1:nbins
    subplot(rows,cols,ibin+cols)
%     subBinRunMean(iSub,iRoi,ibin,rwd,:)
%     subBinRunFft(iSub,iRoi,ibin,rwd,:)
    for rwd=1:2
        plot(1:(trialLength*ntrials)/2, squeeze(binFft(iRoi,ibin,rwd,:)),'.-','Color', plotColors{rwd},'linewidth',linewidth,'markersize',markersize);
        hold on
    end
    xlabel('frequency (cycles/run)');
        t = {'Fourier amplitude'; ['(' fmriUnits ')'] };
        if ibin==1
       ylabel(t);
        end
%     ylabel(['Fourier amplitude  (' fmriUnits ')']);
end
%%
for j=i-1:i
    figure(j)
    subplot(rows,cols,3)
    ylim = get(gca,'ylim');
    yAxisMin = ylim(1);
    yAxisMax = ylim(2);
    yticklabel = get(gca,'yTickLabel');
    ytick = get(gca,'yTick');
    
    for isubplot=1:rows*cols
        subplot(rows,cols,isubplot)
%         set(gca,'ylim',ylim);
%         drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -14/64, 'xAxisMargin', 4/64, 'yAxisMargin', 8/64,'xAxisMinMaxSetByTicks',0,...
%             'labelFontSize',fontsize,'yAxisMin',yAxisMin,'yAxisMax',yAxisMax,'yAxisMinMaxSetByTicks',1,'yTickLabel',yticklabel,'yTick',ytick);
        drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -14/64, 'xAxisMargin', 4/64, 'yAxisMargin', 8/64,'xAxisMinMaxSetByTicks',0,...
            'labelFontSize',fontsize,'yTickLabel',yticklabel,'yTick',ytick);

        %    if isubplot<4
        %        xlabel('time (sec)');
        %        ylabel(['BOLD signal (' fmriUnits ')']);
        %    else
        %        xlabel('frequency (cycles/run)');
        %        ylabel(['Fourier amplitude  (' fmriUnits ')']);
        %    end
        axis square
    end
    set(gcf,'position',[10 10 20.5 	14]);
    
end
figure(i-1); print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig3_3bins_' ConcatProjStr '.pdf']);
figure(i); print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig3_3bins_group_' ConcatProjStr '.pdf']);

%%
i=i+1; figure(i); clf
dx=0.1;
xvec=-15:dx:15;
yvec=-15:dx:15;
[X Y] = meshgrid(xvec,yvec);
[TH R] = cart2pol(X, Y);
R(R>15) = 15;
imagesc(R);
colormap('cool');
xticks([]); yticks([]);
axis image
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig2_circular_cmap.pdf']);




%%
toc

% save([dataFolder 'rwdEcc_concat' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
%     'subFolders', 'roiNames', ...
%     'expName','stairThresh','trialLength',...
%     'globalMean','regressBetasGlobal','runRwd', ...
%     'pVal_rfx', 'pVal_ffx','meanRealDiff','meanRealDiffStd','nperms','binBorders','nbins','numBinVoxels',...
%     'pVal_var_timecourse','pVal_var','pVal_ampVar','pVal_phVar','pVal_otherVar');

