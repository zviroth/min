close all; clear all;
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean=0;

iRoi=2;
rows=1;
cols=3;
subplots = {1,2,3,4};
% rows=2;
% cols = 9;
% subplots = {1:3, 4:6, 7:9, cols+1:cols+3, cols+4:cols+6, cols+7:cols+9};
% subplots = {1 3 5};
linelength = 0.2;
linewidth = 1;
markersize = 10;
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
fontsize=9;
dsSurfaceContrast = 1;
    dsSurfaceAlpha = 0.15;
   
    
    i=0; 
i=i+1; figure(i); clf



saveFolder = '/Volumes/MH02086153MACDT-Drobo/allMinBehavioral/';
onlyCorrect=0;%1=correct,2=incorrect,0=all (with response)
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
end
tic
load([saveFolder 'behavioralData' onlyCorrectString '.mat'], 'subFolders', 'subPupil', 'runSize', 'rwdPupil',...
    'meanPupil','stdRwd','diffPupil','expName',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh','subRTvar','rwdLevel','numTrials',...
    'trialCorrectness' ,'trialResponse','trialRT','propCorrect','stairThresh');
plotColors = { [1 0 0],[0 0 1],[0 1 0], [0.5 1 0.2]};

%%
minLength = 2000;
for rwd=1:2
    groupAll{rwd}(1,:) = meanPupil{1,rwd}(1:minLength);
    for iSub=2:length(subFolders)
        groupAll{rwd}(iSub,:) =  meanPupil{iSub,rwd}(1:minLength);
    end
    groupMean{rwd} = mean(groupAll{rwd});
    groupStd{rwd} = std(groupAll{rwd});
end
clear diffPupil groupAll
for iSub = 1:length(subFolders)
%     minlength = min(length(meanPupil{iSub,1}), length(meanPupil{iSub,2}));
    diffPupil(iSub,:) = meanPupil{iSub,1}(1:minLength) - meanPupil{iSub,2}(1:minLength);
end

% iSub=6;
% subplot(rows,cols,1)
% for rwd=1:2
%     dsErrorsurface((1:minLength)*2, meanPupil{iSub,rwd}(1:minLength), stdRwd{iSub,rwd}(1:minLength)./sqrt(numTrials(iSub,rwd)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
%     hold on
% end
% for rwd=1:2
%     plot((1:minLength)*2,meanPupil{iSub,rwd}(1:minLength),'Color',plotColors{rwd},'linewidth',linewidth);
% end
% subplot(rows,cols,2)
% plot((1:minLength)*2, diffPupil(iSub,:),'k','linewidth',linewidth);
% hline(0);

subplot(rows,cols,subplots{1})
for rwd=1:2
    dsErrorsurface((1:minLength)*2, groupMean{rwd}(1:minLength), groupStd{rwd}(1:minLength)./sqrt(length(subFolders)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
end
for rwd=1:2
    plot((1:minLength)*2,groupMean{rwd}(1:minLength),'Color',plotColors{rwd},'linewidth',linewidth);
end

% subplot(rows,cols,subplots{2})
% meanDiff = mean(diffPupil);
% stdDiff = std(diffPupil);
% dsErrorsurface((1:minLength)*2, meanDiff(1:minLength), stdDiff(1:minLength)./sqrt(length(subFolders)), 'k',dsSurfaceAlpha);
% hold on
% plot((1:minLength)*2, meanDiff(1:minLength),'k','linewidth',linewidth);
% hline(0);




%% HEART RATE

% ConcatProj = 0;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';

onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
elseif onlyCorrect==0
    onlyCorrectString = '_validresponse';
end
zScoreString = '';
fmriUnits = '% change image intensity';
if toZscore
    zScoreString = '_zscored';
    fmriUnits = 'std image intensity';
end

globalMeanString = '';
if regressGlobalMean
    globalMeanString = '_globalRegressed';
end

for ConcatProj=0:1
    clear subPulseTrial numBinVoxels binMeanTserie subBinTrialResponse
    clear binPulseKernel binPulseResidualTC binPhysioKernel binPhysioResidualTC binRespKernel binRespResidualTC
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
        'allGoodTrials','allGoodTRs');
    

    % ROIs = 1:length(roiNames)-1;
    ROIs = 1:length(roiNames);
    % ROIs = [2];
    eccMin = 0.2;
    eccMax = 70;
    nbins = 12;
    TR=1.5;
    binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
    nbins = length(binBorders)-1;
    for ibin=2:length(binBorders)
        binCenters(ibin-1) = (binBorders(ibin)+binBorders(ibin-1))/2;
    end
    
    goodSubs = [1:length(subFolders)];
    
    
    
    for iSub = 1:length(goodSubs)%1:length(goodSubs)%length(subdirs)
        %     iSub = goodSubs(i)
        for iRoi=1:length(roiNames)
            for rwd=1:2
                
                trRV{iSub,rwd} = reshape(rwdRvTC{goodSubs(iSub),rwd}, ecgTrial,[]);%this is only good trials!
                subMeanRV(iSub,rwd,:) = nanmean(trRV{iSub,rwd},2);
                subRespStd(iSub,rwd) = std(subMeanRV(iSub,rwd,:));%std amplitude of mean
                subRespVar(iSub,rwd,:) = nanstd(trRV{iSub,rwd},0,2);%timepoint variability
                respStd = nanstd(trRV{iSub,rwd});%std amp per trial
                subRespStdVar(iSub,rwd) = std(respStd);
                f = fft(trRV{iSub,rwd});
                subRespPhVar(iSub,rwd) = nanstd(angle(f(2,:)));
                subRespAmpVar(iSub,rwd) = nanstd(abs(f(2,:)));
                
                
                trPulse{iSub,rwd} = reshape(rwdPulseTC{goodSubs(iSub),rwd}, ecgTrial,[]);%this is only good trials!
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
                designMatPulse{goodSubs(iSub),rwd}= designMatPulse{goodSubs(iSub),rwd}*ecgSampleRate;%change to beats per sec
                pulseKernel(iSub,iRoi,rwd,:) = designMatPulse{goodSubs(iSub),rwd}'\subTrialResponse{goodSubs(iSub),iRoi,rwd}(:);%this is only good trials!
                pulseResidualTC{iSub,iRoi,rwd} = subTrialResponse{goodSubs(iSub),iRoi,rwd}(:)' - squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{goodSubs(iSub),rwd};
                
                designMatResp{goodSubs(iSub),rwd} = designMatResp{goodSubs(iSub),rwd}*ecgSampleRate;%change to beats per sec
                respKernel(iSub,iRoi,rwd,:) = designMatResp{goodSubs(iSub),rwd}'\subTrialResponse{goodSubs(iSub),iRoi,rwd}(:);
                respResidualTC{iSub,iRoi,rwd} = subTrialResponse{goodSubs(iSub),iRoi,rwd}(:)' - squeeze(respKernel(iSub,iRoi,rwd,:))'*designMatResp{goodSubs(iSub),rwd};
                
                designMatRespPulse{goodSubs(iSub),rwd} = designMatRespPulse{goodSubs(iSub),rwd}*ecgSampleRate;%change to beats per sec
                physioKernel(iSub,iRoi,rwd,:) = designMatRespPulse{goodSubs(iSub),rwd}'\subTrialResponse{goodSubs(iSub),iRoi,rwd}(:);
                physioResidualTC{iSub,iRoi,rwd} = subTrialResponse{goodSubs(iSub),iRoi,rwd}(:)' - squeeze(physioKernel(iSub,iRoi,rwd,:))'*designMatRespPulse{goodSubs(iSub),rwd};
                subTrialResponse{goodSubs(iSub),iRoi,rwd} = reshape(physioResidualTC{iSub,iRoi,rwd}(:), trialLength, length(physioResidualTC{iSub,iRoi,rwd}(:))/trialLength);
                
                
            end
        end
    end
    
    % regress physio per bin

%     for iSub = 1:length(goodSubs)%1:length(goodSubs)%length(subdirs)
%         for rwd=1:2
%             for ibin=1:nbins
%                 binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
%                 numBinVoxels(iSub,ibin) = sum(binVoxels);
%                 %                 binVoxels = binVoxels & areas{iSub,iRoi}==1;%ONLY V1
%                 binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
%                 binMeanTseries{iSub,ibin,rwd} = nanmean(roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,allGoodTRs{iSub,rwd}));%mean timecourse across voxels
%                 subBinTrialResponse{iSub,ibin,rwd} = reshape(binMeanTseries{iSub,ibin,rwd}, trialLength, length(binMeanTseries{goodSubs(iSub),ibin,rwd})/trialLength);
%                 
%                 binPulseKernel(iSub,ibin,rwd,:) = designMatPulse{goodSubs(iSub),rwd}'\binMeanTseries{iSub,ibin,rwd}';
%                 binPulseResidualTC{iSub,ibin,rwd} = binMeanTseries{iSub,ibin,rwd} - squeeze(binPulseKernel(iSub,ibin,rwd,:))'*designMatPulse{goodSubs(iSub),rwd};
%                 
%                 binPhysioKernel(iSub,ibin,rwd,:) = designMatRespPulse{goodSubs(iSub),rwd}'\binMeanTseries{iSub,ibin,rwd}';
%                 binPhysioResidualTC{iSub,ibin,rwd} = binMeanTseries{iSub,ibin,rwd} - squeeze(binPhysioKernel(iSub,ibin,rwd,:))'*designMatRespPulse{goodSubs(iSub),rwd};
%                 
%                 binRespKernel(iSub,ibin,rwd,:) = designMatResp{goodSubs(iSub),rwd}'\binMeanTseries{iSub,ibin,rwd}';
%                 binRespResidualTC{iSub,ibin,rwd} = binMeanTseries{iSub,ibin,rwd} - squeeze(binRespKernel(iSub,ibin,rwd,:))'*designMatResp{goodSubs(iSub),rwd};
%             end
%         end
%     end
    iRoi = 2;
    
%     realSubPulseKernelStd = std(pulseKernel,0,4);
%     realSubRespKernelStd = std(respKernel,0,4);
%     pulseKernelStd = std(binPulseKernel,0,4);
%     meanPulseKernelStd = squeeze(mean(pulseKernelStd));
%     stdPulseKernelStd = squeeze(std(pulseKernelStd));
%     pulseKernelStdDiff = squeeze(pulseKernelStd(:,:,1) - pulseKernelStd(:,:,2));
%     meanPulseKernelStdDiff = squeeze(mean(pulseKernelStdDiff));
%     stdPulseKernelStdDiff = squeeze(std(pulseKernelStdDiff));
    
    %heart rate
    for iSub=1:length(subFolders)
        for rwd=1:2
            temp = reshape(rwdPulseTC{iSub,rwd},ecgTrial,[]);
            subPulseTrial(iSub,rwd,:) = nanmean(temp,2);
        end
    end
    subPulseTrial = subPulseTrial*ecgSampleRate;%changing to beats-per-second
    meanPulseTrial = squeeze(nanmean(subPulseTrial));
    stdPulseTrial = squeeze(nanstd(subPulseTrial));
    subPulseDiff = squeeze(subPulseTrial(:,1,:) - subPulseTrial(:,2,:));
    meanPulseDiff = squeeze(nanmean(subPulseDiff));
    stdPulseDiff = squeeze(nanstd(subPulseDiff));
    
    subplot(rows,cols,subplots{3})
    if ConcatProj==0
        lineColor = [0 0 0];
        linestyle = '-';
    else
        lineColor = [0 0 0];
        linestyle = '--';
    end
    meanRwdPulseKernel = squeeze(mean(pulseKernel,3));%average across rwd, per subject
    meanPulseKernel = squeeze(mean(meanRwdPulseKernel(:,iRoi,:)));%average across subjects
    stdPulseKernel = squeeze(std(meanRwdPulseKernel(:,iRoi,:)));
    dsErrorsurface(TR*(1:trialLength), squeeze(meanPulseKernel), squeeze(stdPulseKernel)./sqrt(size(pulseKernel,1)), dsSurfaceContrast*lineColor,dsSurfaceAlpha);
    hold on
    plot(TR*(1:trialLength), squeeze(meanPulseKernel),linestyle,'Color', lineColor,'linewidth',linewidth,'markersize',markersize);
    
    
%     t = {'BOLD signal'; ['(' fmriUnits ')'] };
%     ylabel(t);
%     ylim([-.02 0.03]);
    % ylabel('pulse kernel');
    %     xlabel('time (sec)');
    %     drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -14/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',1,...
    %         'labelFontSize',fontsize);
    
    %KERNEL AMPLITUDE
    kernelAmp(ConcatProj+1) = std(squeeze(meanPulseKernel));
    kernelAmpSub(ConcatProj+1,:) = std(squeeze(meanRwdPulseKernel(:,iRoi,:))');
    
end
t = {'BOLD signal'; ['(' fmriUnits ')'] };
    ylabel(t);
kernelAmp(2)/kernelAmp(1);
kernelAmpSub(2,:)'./kernelAmpSub(1,:)';

% mean pulse rate
subplot(rows,cols,subplots{2})
for rwd=1:2
    dsErrorsurface((0:ecgTrial-1)/ecgSampleRate, meanPulseTrial(rwd,:),stdPulseTrial(rwd,:)./sqrt(size(rwdPulseTC,1)), dsSurfaceContrast*plotColors{rwd},dsSurfaceAlpha);
    hold on
    plot((0:ecgTrial-1)/ecgSampleRate,meanPulseTrial(rwd,:), 'Color', plotColors{rwd}, 'linewidth', linewidth,'markersize',20);
end
ylabel('pulse (beats/sec)');





% %% BOLD RESPONSE
% subplot(rows,cols,subplots{4})
% for ConcatProj = 0:1
%     if ConcatProj==0
%         lineColor = [0 0 0];
%         linestyle = '-';
%     else
%         lineColor = [0 0 0];
%         linestyle = '--';
%     end
%     ConcatProjStr = '';
%     if ConcatProj
%         ConcatProjStr = 'ConcatProj';
%     end
%     load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'subResponse');
%     subRoiResponse{ConcatProj+1} = squeeze(mean(subResponse(:,2,:,:),3));
%     meanResponse = squeeze(mean(subRoiResponse{ConcatProj+1}));    
%     stdResponse = squeeze(std(subRoiResponse{ConcatProj+1}));
%     hold on
%     dsErrorsurface(TR*(0:trialLength-1), meanResponse, stdResponse./sqrt(size(subResponse,1)), dsSurfaceContrast*lineColor,dsSurfaceAlpha);
%     hold on
%     plot(TR*(0:trialLength-1), meanResponse,linestyle,'Color', lineColor,'linewidth',linewidth,'markersize',markersize);
% end
% t = {'BOLD signal'; ['(' fmriUnits ')'] };
%     ylabel(t);
%%

for isubplot=1%:2
   subplot(rows,cols,subplots{isubplot})
   xlabel('time (ms)');
   axis square
   if mod(isubplot,2)>0
       ylabel('pupil size (arb. units)');
        drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -14/64, 'xAxisMargin', 4/64, 'yAxisMargin', 4/64,'labelFontSize',fontsize);
   else
       ylabel('\Delta pupil size (arb. units)');
       drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -12/64, 'xAxisMargin', 4/64, 'yAxisMargin', 2/64,'labelFontSize',fontsize);
   end
end
for isubplot=2:3
    subplot(rows,cols,subplots{isubplot})
    xlabel('time (s)');
    yticklabel = get(gca,'yTickLabel');
    ytick = get(gca,'yTick');
    if isubplot==3
        yticklabel = {'-0.015','0','0.03'};
        ytick = [-0.015 0 0.03];
%         yticklabel = yticklabel(1:2:end);
%         ytick = ytick(1:2:end);
    end
    drawPublishAxis('xLabelOffset', -8/64,'yLabelOffset', -16/64, 'xAxisMargin', 8/64, 'yAxisMargin', 4/64,'xAxisMinMaxSetByTicks',1,...
        'labelFontSize',fontsize,'yTickLabel',yticklabel,'yTick',ytick);
    axis square
end

% get(gca,'yTickLabel')


%% SAVE
% set(gcf,'position',[10 10 21 18]);
set(gcf,'position',[10 10 21 8]);
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig2_PupilAndPulse.pdf']);
%
