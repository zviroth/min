%first run associatePhysio_min.m on the not-concatenated folders
%then copy all stimfiles into the concatenated folders.
% dataFolder = '~/data/rwdFmri/';
close all
mrQuit
clear all
tic
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
% dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects/';
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', ...
%     '003220180105', '0034i20180209', '003520180328', '0039i20180220', '004020180328','004120180320', ...
%     '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
subFolders = {'000520180116', '0008i20180213', '0016i20180207',  ...
    '003220180105', '0034i20180209', '0039i20180220', '004020180328','004120180320', ...
    '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
% subFolders = {'000520180116',	'004120180320','0008i20180213','0042i20180412',...
%     'extra_000820171219','0016i20180207','0039i20180220','0045i20180309','extra_001620171213',...
%     '002220171212','004020180221','0046i20180409','003220180105','004020180328','0049i20180404',...
%     '0034i20180209','004120180223','005220180621'};

roiNames = {'leftBenson','rightBenson','lh_17Networks_16','rh_17Networks_16','leftCerebellarCortex','rightCerebellarCortex'};

cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end
junkedFrames = 10;
trialLength=10;
clear ecg resp
TR=1.5;
ecgselect=0.009;
respselect = 0.1;
display=0;
ecgSampleRate = 50;
ecgTrial = ecgSampleRate*TR*trialLength;
trialsPerRun=15;
ecgRunLength = ecgTrial*(trialsPerRun+1);
ecgInterpMethod = 'linear';
respInterpMethod = 'linear';
respStdWindow = 6*ecgSampleRate;



deconvLength = 12;
numSubs = length(subFolders);
for iSub=1:numSubs
    cd(subFolders{iSub});
    v=newView;
% 
    concatGroupNum = viewGet(v,'groupNum','concatSessions'); %concatenation of multiple sessions
    if isempty(concatGroupNum)%single session
        concatGroupNum = viewGet(v,'groupNum','Concatenation');
    end
%     concatGroupNum = viewGet(v,'groupNum','Concatenation');
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', concatGroupNum);

    nScans = viewGet(v, 'nscans');
    clear params

    for iscan = 1:2%nScans%2 concatenations, 1 for each reward type
        s = viewGet(v, 'stimfile', iscan);
        rwdType = s{1}.stimulus.rewardVal;
        if strcmp(rwdType, 'H')
            rwd = 1;
        elseif strcmp(rwdType, 'L')
            rwd = 2;
        else
            keyboard
        end
        numRuns(iSub,rwd) = length(s);
        concatInfo{iSub,rwd} = viewGet(v, 'concatInfo', iscan);
        numTRs(iSub,rwd) = length(concatInfo{iSub,rwd}.whichScan);
        
        %get the run #s for this concatenation. Then locate their ECG
        %traces. 

        %load all voxels timecourse
        origTSeries = loadTSeries(v,iscan);
        origSize{iSub,rwd} = size(origTSeries);
        allTseries{iSub,rwd} = reshape(origTSeries,[],numTRs(iSub,rwd));%voxels,TRs
%         nullTseries{iSub,rwd} = allTseries{iSub,rwd}(:,nullTrialsTRs==1);
%         stimTseries{iSub,rwd} = allTseries{iSub,rwd}(:,nullTrialsTRs==0);
        globalMean{iSub,rwd} = nanmean(allTseries{iSub,rwd})';

        
        for r=1:numRuns(iSub,rwd)
            
            %RESPIRATION
            respfilename = fullfile(s{r}.myscreen.originalmlrdir,s{r}.myscreen.physiodir,s{r}.myscreen.respirationfilename);
            if ~isfile(respfilename)
                %if the directory has been moved, but this is the original
                %session this should work.
                respfilename = fullfile(s{r}.myscreen.physiodir,s{r}.myscreen.respirationfilename);
                if ~isfile(respfilename)%concatenated session, and original directory has been moved?
                    keyboard
                end
            end
            resp{iSub,rwd,r}=load(respfilename);
            [respPeaks{iSub,rwd,r},criterion] = pickpeaks(resp{iSub,rwd,r},respselect,display);
            respPeaksDiff{iSub,rwd,r} = diff(respPeaks{iSub,rwd,r});
            respPeaksAmp{iSub,rwd,r} = resp{iSub,rwd,r}(respPeaks{iSub,rwd,r});
            interpPeaksAmp{iSub,rwd,r} = interp1(respPeaks{iSub,rwd,r}, respPeaksAmp{iSub,rwd,r}, 1:length(resp{iSub,rwd,r}), respInterpMethod, 'extrap');
            interpPeaksAmp{iSub,rwd,r}(end:ecgRunLength) = NaN;
            rv{iSub,rwd,r} = zeros(size(resp{iSub,rwd,r}));
            
            [respTroughs{iSub,rwd,r},criterion] = pickpeaks(-resp{iSub,rwd,r},respselect,display);%gives the times
            respTroughsAmp{iSub,rwd,r} = resp{iSub,rwd,r}(respTroughs{iSub,rwd,r});
            interpTroughsAmp{iSub,rwd,r} = interp1(respTroughs{iSub,rwd,r}, respTroughsAmp{iSub,rwd,r}, 1:length(resp{iSub,rwd,r}), respInterpMethod, 'extrap');
            interpTroughsAmp{iSub,rwd,r}(end:ecgRunLength) = NaN;

            respRateTime{iSub,rwd,r} = respPeaks{iSub,rwd,r}(1:end-1) + 0.5*diff(respPeaks{iSub,rwd,r});%timepoints between the peaks
            interpRespPeaksDiff{iSub,rwd,r} = interp1(respRateTime{iSub,rwd,r}, respPeaksDiff{iSub,rwd,r}, 1:length(resp{iSub,rwd,r}), respInterpMethod, 'extrap');
            interpRespPeaksDiff{iSub,rwd,r}(end:ecgRunLength) = NaN;
            
            rvt{iSub,rwd,r} = (interpPeaksAmp{iSub,rwd,r} - interpTroughsAmp{iSub,rwd,r})./interpRespPeaksDiff{iSub,rwd,r};
            
            %we will junk the first 10 TRs of BOLD anyway, so starting from
            %6 sec. BUT, this measure is lagging.
            for t=respStdWindow:length(resp{iSub,rwd,r})
                rv{iSub,rwd,r}(t) = nanstd(resp{iSub,rwd,r}(t-respStdWindow+1:t));
            end
            rv{iSub,rwd,r}(1:respStdWindow) = rv{iSub,rwd,r}(t);
            rv{iSub,rwd,r}(end:ecgRunLength) = NaN;
            rvt{iSub,rwd,r}(end:ecgRunLength) = NaN;
            temp = reshape(rv{iSub,rwd,r}, ecgTrial,[]);
%             rv{iSub,rwd,r}(1:respStdWindow) = squeeze(mean(temp(1:respStdWindow,:),2));%extrapolate to first trial from mean across trials
            trialsRV{iSub,rwd}(r,:,:) = reshape(rv{iSub,rwd,r}, ecgTrial,[]);
            trialsRVT{iSub,rwd}(r,:,:) = reshape(rvt{iSub,rwd,r}, ecgTrial,[]);
            meanTrialRvRun{iSub,rwd}(r,:) = squeeze(mean(trialsRV{iSub,rwd}(r,:,2:end),3));%mean across trials, exclude first trial
            meanTrialRvtRun{iSub,rwd}(r,:) = squeeze(mean(trialsRVT{iSub,rwd}(r,:,2:end),3));%mean across trials, exclude first trial
            
            
            
            %ECG
            %first load the ECG file
            ecgfilename = fullfile(s{r}.myscreen.originalmlrdir,s{r}.myscreen.physiodir,s{r}.myscreen.ecgfilename);
            if ~isfile(ecgfilename)
                %if the directory has been moved, but this is the original
                %session this should work.
                ecgfilename = fullfile(s{r}.myscreen.physiodir,s{r}.myscreen.ecgfilename);
                if ~isfile(ecgfilename)%concatenated session, and original directory has been moved?
                    keyboard
                end
            end
            ecg{iSub,rwd,r}=load(ecgfilename);
            [ecgPeaks{iSub,rwd,r},criterion] = pickpeaks(ecg{iSub,rwd,r},ecgselect,display);
            ecgPeaksDiff{iSub,rwd,r} = diff(ecgPeaks{iSub,rwd,r});
            ecgPeaksAmp{iSub,rwd,r} = ecg{iSub,rwd,r}(ecgPeaks{iSub,rwd,r});
            ecgRateTime{iSub,rwd,r} = ecgPeaks{iSub,rwd,r}(1:end-1) + 0.5*diff(ecgPeaks{iSub,rwd,r});%timepoints between the peaks
            scaleFactor = mean(ecgPeaksAmp{iSub,rwd,r})/mean(ecgPeaksDiff{iSub,rwd,r});%for visualization purposes
            interpPeaksDiff{iSub,rwd,r} = interp1(ecgRateTime{iSub,rwd,r}, ecgPeaksDiff{iSub,rwd,r}, 1:length(ecg{iSub,rwd,r}), ecgInterpMethod, 'extrap');
            interpPeaksDiff{iSub,rwd,r}(end:ecgRunLength) = NaN;
            
            trialsPeakDiff{iSub,rwd}(r,:,:) = reshape(interpPeaksDiff{iSub,rwd,r}, ecgTrial,[]);
            meanTrialPeaksDiffRun{iSub,rwd}(r,:) = squeeze(mean(trialsPeakDiff{iSub,rwd}(r,:,:),3));%mean across trials
        
            ecgPulseRate{iSub,rwd,r} = 1./ecgPeaksDiff{iSub,rwd,r};
            interpPulseRate{iSub,rwd,r} = interp1(ecgRateTime{iSub,rwd,r}, ecgPulseRate{iSub,rwd,r}, 1:length(ecg{iSub,rwd,r}), ecgInterpMethod, 'extrap');
            interpPulseRate{iSub,rwd,r}(end:ecgRunLength) = NaN;
            
            trialsPulseRate{iSub,rwd}(r,:,:) = reshape(interpPulseRate{iSub,rwd,r}, ecgTrial,[]);
            meanTrialPulseRateRun{iSub,rwd}(r,:) = squeeze(mean(trialsPulseRate{iSub,rwd}(r,:,:),3));%mean across trials
           
        end
        %RESPIRATION
        rwdMeanRVT(iSub,rwd,:) = mean(meanTrialRvtRun{iSub,rwd});
        rwdMeanRV(iSub,rwd,:) = mean(meanTrialRvRun{iSub,rwd});
        rwdRvTC{iSub,rwd} = [];
        rwdRvtTC{iSub,rwd} = [];
        for r=1:numRuns(iSub,rwd)
            %removing junked frames & concatenate
            rwdRvTC{iSub,rwd} = horzcat(rwdRvTC{iSub,rwd}, rv{iSub,rwd,r}(ecgSampleRate*1.5*junkedFrames+1:end)');
            rwdRvtTC{iSub,rwd} = horzcat(rwdRvtTC{iSub,rwd}, rvt{iSub,rwd,r}(ecgSampleRate*1.5*junkedFrames+1:end));
        end
        %downsample pulse using median
        trRV = reshape(rwdRvTC{iSub,rwd}, ecgSampleRate*1.5,[]);
        trRVT = reshape(rwdRvtTC{iSub,rwd}, ecgSampleRate*1.5,[]);
        downsampledRV{iSub,rwd} = nanmedian(trRV);
        downsampledRVT{iSub,rwd} = nanmedian(trRVT);
        %create design matrix for deconvolution
        designMatRv{iSub,rwd}(1,:) = downsampledRV{iSub,rwd};
        designMatRvt{iSub,rwd}(1,:) = downsampledRVT{iSub,rwd};
        for i=2:deconvLength
            designMatRv{iSub,rwd}(i,:) = circshift(designMatRv{iSub,rwd}(i-1,:),1);
            designMatRvt{iSub,rwd}(i,:) = circshift(designMatRvt{iSub,rwd}(i-1,:),1);
        end
%         designMatResp{iSub,rwd} = [designMatRv{iSub,rwd}; designMatRvt{iSub,rwd}];
        designMatResp{iSub,rwd} = [designMatRv{iSub,rwd}];%ONLY RV
        
        %ECG
        rwdMeanPeaksDiff(iSub,rwd,:) = mean(meanTrialPeaksDiffRun{iSub,rwd});
        rwdMeanPulseRate(iSub,rwd,:) = mean(meanTrialPulseRateRun{iSub,rwd});
        %concatenate all pulse traces for this rwd
        rwdPulseTC{iSub,rwd} = [];
        for r=1:numRuns(iSub,rwd)
            %removing junked frames & concatenate
            rwdPulseTC{iSub,rwd} = horzcat(rwdPulseTC{iSub,rwd}, interpPulseRate{iSub,rwd,r}(ecgSampleRate*1.5*junkedFrames+1:end));
        end
        %downsample pulse using median
        trPulse = reshape(rwdPulseTC{iSub,rwd}, ecgSampleRate*1.5,[]);
        downsampledPulse{iSub,rwd} = nanmedian(trPulse);
        %create design matrix for deconvolution
        designMatPulse{iSub,rwd}(1,:) = downsampledPulse{iSub,rwd};
        for i=2:deconvLength
            designMatPulse{iSub,rwd}(i,:) = circshift(designMatPulse{iSub,rwd}(i-1,:),1);
        end
        
        designMatRespPulse{iSub,rwd} = [designMatResp{iSub,rwd}; designMatPulse{iSub,rwd}];
        
        %load fMRI timeseries
        for iRoi = 1:length(roiNames)
           roiTC{iSub,iRoi,rwd} = loadROITSeries(v, roiNames{iRoi}, iscan, [], 'keepNAN',true);
           roiTC{iSub,iRoi,rwd}.tSeries = zscoreConcat(roiTC{iSub,iRoi,rwd}.tSeries, concatInfo{iSub,rwd});
           
           regressBetasGlobal{iSub,rwd,iRoi} = globalMean{iSub,rwd}\roiTC{iSub,iRoi,rwd}.tSeries';
           roiTCregress{iSub,iRoi,rwd}.tSeries = roiTC{iSub,iRoi,rwd}.tSeries - regressBetasGlobal{iSub,rwd,iRoi}'*globalMean{iSub,rwd}';
           
        end
        
        
        for iRoi = 1:length(roiNames)
           meanRoiTC{iSub,iRoi,rwd} = nanmean(roiTC{iSub,iRoi,rwd}.tSeries);
           meanRoiTCregress{iSub,iRoi,rwd} = nanmean(roiTCregress{iSub,iRoi,rwd}.tSeries);
           origVar(iSub,iRoi,rwd) = var(meanRoiTC{iSub,iRoi,rwd});
           origVarRegressGlobal(iSub,iRoi,rwd) = var(meanRoiTCregress{iSub,iRoi,rwd});

           pulseBoldCorr(iSub,iRoi,rwd) = corr(downsampledPulse{iSub,rwd}',meanRoiTC{iSub,iRoi,rwd}');
           
           pulseKernel(iSub,iRoi,rwd,:) = designMatPulse{iSub,rwd}'\meanRoiTC{iSub,iRoi,rwd}';
           pulseKernelRegress(iSub,iRoi,rwd,:) = designMatPulse{iSub,rwd}'\meanRoiTCregress{iSub,iRoi,rwd}';
           pulseResidualTC{iSub,iRoi,rwd} = meanRoiTC{iSub,iRoi,rwd} - squeeze(pulseKernel(iSub,iRoi,rwd,:))'*designMatPulse{iSub,rwd};
           pulseResidualVar(iSub,iRoi,rwd) = var(pulseResidualTC{iSub,iRoi,rwd});
           pulseResidualTCregressGlobal{iSub,iRoi,rwd} = meanRoiTCregress{iSub,iRoi,rwd} - squeeze(pulseKernelRegress(iSub,iRoi,rwd,:))'*designMatPulse{iSub,rwd};
           pulseResidualVarRegressGlobal(iSub,iRoi,rwd) = var(pulseResidualTCregressGlobal{iSub,iRoi,rwd});
           
           physioKernel(iSub,iRoi,rwd,:) = designMatRespPulse{iSub,rwd}'\meanRoiTC{iSub,iRoi,rwd}';
           physioKernelRegress(iSub,iRoi,rwd,:) = designMatRespPulse{iSub,rwd}'\meanRoiTCregress{iSub,iRoi,rwd}';
           physioResidualTC{iSub,iRoi,rwd} = meanRoiTC{iSub,iRoi,rwd} - squeeze(physioKernel(iSub,iRoi,rwd,:))'*designMatRespPulse{iSub,rwd};
           physioResidualVar(iSub,iRoi,rwd) = var(physioResidualTC{iSub,iRoi,rwd});
           physioResidualTCregressGlobal{iSub,iRoi,rwd} = meanRoiTCregress{iSub,iRoi,rwd} - squeeze(physioKernelRegress(iSub,iRoi,rwd,:))'*designMatRespPulse{iSub,rwd};
           physioResidualVarRegressGlobal(iSub,iRoi,rwd) = var(physioResidualTCregressGlobal{iSub,iRoi,rwd});
           
           respKernel(iSub,iRoi,rwd,:) = designMatResp{iSub,rwd}'\meanRoiTC{iSub,iRoi,rwd}';
           respKernelRegress(iSub,iRoi,rwd,:) = designMatResp{iSub,rwd}'\meanRoiTCregress{iSub,iRoi,rwd}';
           respResidualTC{iSub,iRoi,rwd} = meanRoiTC{iSub,iRoi,rwd} - squeeze(respKernel(iSub,iRoi,rwd,:))'*designMatResp{iSub,rwd};
           respResidualVar(iSub,iRoi,rwd) = var(respResidualTC{iSub,iRoi,rwd});
           respResidualTCregressGlobal{iSub,iRoi,rwd} = meanRoiTCregress{iSub,iRoi,rwd} - squeeze(respKernelRegress(iSub,iRoi,rwd,:))'*designMatResp{iSub,rwd};
           respResidualVarRegressGlobal(iSub,iRoi,rwd) = var(respResidualTCregressGlobal{iSub,iRoi,rwd});
           
        end
    end

    deleteView(v);
    
    
    %% return to home directory
    cd('..');
end

pulseRsqr = 1 - pulseResidualVar./origVar;
pulseRsqrRegressGlobal = 1 - pulseResidualVarRegressGlobal./origVarRegressGlobal;
globalMeanRsqr = 1 - origVarRegressGlobal./origVar;
%%
save([dataFolder 'rwd_physioRegress.mat'], 'dataFolder', 'subFolders', 'numSubs','roiNames', ...
    'numRuns','numTRs','concatInfo',...
    'trialsPerRun','trialLength','junkedFrames','TR',...
    'ecgselect','ecgSampleRate','ecgTrial','ecgRunLength','ecgInterpMethod',...
    'numRuns','concatInfo','numTRs','globalMean',...
    'ecg','ecgPeaksDiff','ecgPeaksAmp','ecgRateTime','interpPeaksDiff',...
    'trialsPeakDiff','ecgPulseRate','interpPulseRate',...
    'trialsPulseRate','meanTrialPulseRateRun','meanTrialPeaksDiffRun',...
    'rwdMeanPeaksDiff','rwdMeanPulseRate','rwdPulseTC','downsampledPulse',...
    'rwdRvTC','rwdRvtTC','downsampledRV','downsampledRVT',...
    'origVar','origVarRegressGlobal','pulseResidualVar','pulseResidualVarRegressGlobal',...
    'meanRoiTC','meanRoiTCregress',...
    'pulseRsqr','pulseRsqrRegressGlobal',...
    'pulseKernel','pulseKernelRegress',...
    'pulseResidualTC', 'pulseResidualTCregressGlobal',...
    'rv','rvt','respselect','resp','respPeaks','respPeaksDiff',...
    'respPeaksAmp','interpPeaksAmp','respRateTime',...
    'interpRespPeaksDiff','respTroughs','respTroughsAmp','interpTroughsAmp',...
    'physioKernel','physioKernelRegress','physioResidualTC','physioResidualVar',...
    'physioResidualTCregressGlobal','physioResidualVarRegressGlobal',...
    'designMatRespPulse','designMatResp','designMatPulse','rwdMeanRVT','rwdMeanRV',...
    'pulseBoldCorr','deconvLength',...
    'respKernel','respKernelRegress','respResidualTC','respResidualVar',...
    'respResidualTCregressGlobal','respResidualVarRegressGlobal',...
    'downsampledRV','downsampledRVT');
% save([dataFolder 'rwdRapidData_nullTrials.mat'], 'dataFolder', 'subFolders', 'roiNames', ...
%     'numRuns','numTRs','concatInfo',...
%     'frames', 'junkedFrames', 'TR', 'trialsPerRun', 'trialLength', 'nVolumes',...
%     'nullTseries','stimTseries','nullTseriesRegress','stimTseriesRegress',...
%     'rwdTrialPulse','meanRwdTrialPulse');
toc




%%
function newConcat = zscoreConcat(concatData, concatInfo)%voxels,TRs
newConcat = [];
for iRun=1:size(concatInfo.runTransition,1)
    thisRun = concatData(:,concatInfo.runTransition(iRun,1):concatInfo.runTransition(iRun,2));

    newConcat = [newConcat zscore(thisRun,0,2)];
end

end



