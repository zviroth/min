close all
clear all
tic
maxRT=4000;
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean = 0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', ...
%     '003220180105', '0034i20180209', '003520180328', '0039i20180220', '004020180328','004120180320', ...
%     '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
%without s22 and s35
subFolders = {'000520180116', '0008i20180213', '0016i20180207', ...
    '003220180105', '0034i20180209',  '0039i20180220', '004020180328','004120180320', ...
    '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
% subFolders = {'000520180116', '0008i20180213'};
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', ...
%     '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', ...
%     '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
numSubs=length(subFolders);
% roiNames = {'rV1_eccen8', 'lV1_eccen8','rV2_eccen8','lV2_eccen8','rV3_eccen8','lV3_eccen8','rh_17Networks_16','lh_17Networks_16'};
% roiNames = {'leftBenson', 'rightBenson'};
% roiNames = { 'leftBenson', 'rightBenson','lV1_eccen8','rV1_eccen8','lV2_eccen8','rV2_eccen8','lV3_eccen8','rV3_eccen8','lh_17Networks_16','rh_17Networks_16'};
roiNames = { 'leftBenson', 'rightBenson','lh_17Networks_16','rh_17Networks_16'};

junkedFrames = 10;
trialLength=10;
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
deconvLength = 10;
cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end

clear d concatInfo sumResponse subResponse roiMeanTseries meanResponse roiTC subRunResponse trialCorrectness trialResponse trialRT propCorrect subMeanRunTC
regressBetasGlobal={};
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2],[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.'};
rows=2;
cols=ceil(numSubs/rows);
for iSub = 1:numSubs
    
    cd(subFolders{iSub});
    v=newView;
    
    % Find the correct Concatenation group
    if ConcatProj
        concatGroupNum = viewGet(v,'groupNum','ConcatenationProj'); %concatenation of multiple sessions
    else
        
        concatGroupNum = viewGet(v,'groupNum','concat'); %concatenation of multiple sessions
        if isempty(concatGroupNum)%single session
            concatGroupNum = viewGet(v,'groupNum','Concatenation');
        end
    end
    
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', concatGroupNum);
    nScans = viewGet(v, 'nscans');
    clear rois
    
    for iScan = 1:2%nScans%2 concatenations, 1 for each reward type
        s = viewGet(v, 'stimfile', iScan);
        if ~isfield(s{1}, 'stimulus')
            s=s{1};
        end
        rwdType = s{1}.stimulus.rewardVal;
        if strcmp(rwdType, 'H')
            rwd = 1;
        elseif strcmp(rwdType, 'L')
            rwd = 2;
        else
            disp('wtf');
            keyboard
        end
        numRuns(iSub,rwd) = length(s);
        concatInfo{iSub,rwd} = viewGet(v, 'concatInfo', iScan);
        
        %get global mean
        ts = loadTSeries(v,iScan);
        allTseries = reshape(ts,[],size(ts,4));%voxels,TRs
        if toZscore
            allTseries = zscoreConcat(allTseries, concatInfo{iSub,rwd});
        end
        globalMean{iSub,rwd} = nanmean(allTseries)';%(vox,T)
        
        
        
        %BEHAVIOR
        for r=1:length(s)
            %if subject stopped responding at a certain point we
            %need to fill in the last incorrect trials:
            runRwd{iSub,iScan}(r) = s{r}.stimulus.rewardVal;
            trialCorrectness{iSub,rwd}(r,:) = [s{r}.stimulus.trialCorrectness zeros(1,17-size(s{r}.stimulus.trialCorrectness,2))];
            trialResponse{iSub,rwd}(r,:) = [s{r}.stimulus.trialResponse zeros(1,17-size(s{r}.stimulus.trialResponse,2))];%0 means no response
            trialRT{iSub,rwd}(r,:) = [s{r}.stimulus.trialRT NaN(1,17-size(s{r}.stimulus.trialRT,2))];
            propCorrect{iSub,rwd}(r) = s{r}.stimulus.percentCorrect;
            stairThresh{iSub,rwd}(r,:) = [s{r}.stimulus.stair{1}.threshold s{r}.stimulus.stair{2}.threshold];
        end
        subMeanCorrectness(iSub,rwd) = mean(trialCorrectness{iSub,rwd}(:));
        subMeanRT(iSub,rwd) = mean(trialRT{iSub,rwd}(trialRT{iSub,rwd}(:)>0));
        subMedianRT(iSub,rwd) = median(trialRT{iSub,rwd}(:));
        subMeanThresh(iSub,rwd) = mean(stairThresh{iSub,rwd}(:));
        expName{iSub,rwd} = s{iScan}.task{1}.taskFilename;

        temp = trialCorrectness{iSub,rwd}(:,2:end-1);
        trialCorrectnessVec = temp(:);
        temp = trialResponse{iSub,rwd}(:,2:end-1);
        trialResponseVec = temp(:);
        temp = trialRT{iSub,rwd}(:,2:end-1);
        trialRTvec = temp(:);
            
        %GOOD TRIALS
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
        %             trialRT{iSub,rwd};

        goodTRs = repmat(goodTrials',trialLength,1);
        goodTRs = goodTRs(:);
        allGoodTRs{iSub,rwd} = goodTRs;
        allGoodTrials{iSub,rwd} = goodTrials;
        
        
        %get physio data
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
            
            %             [respTroughs{iSub,rwd,r},criterion] = pickpeaks(-resp{iSub,rwd,r},respselect,display);%gives the times
            %             respTroughsAmp{iSub,rwd,r} = resp{iSub,rwd,r}(respTroughs{iSub,rwd,r});
            %             interpTroughsAmp{iSub,rwd,r} = interp1(respTroughs{iSub,rwd,r}, respTroughsAmp{iSub,rwd,r}, 1:length(resp{iSub,rwd,r}), respInterpMethod, 'extrap');
            %             interpTroughsAmp{iSub,rwd,r}(end:ecgRunLength) = NaN;
            %
            %             respRateTime{iSub,rwd,r} = respPeaks{iSub,rwd,r}(1:end-1) + 0.5*diff(respPeaks{iSub,rwd,r});%timepoints between the peaks
            %             interpRespPeaksDiff{iSub,rwd,r} = interp1(respRateTime{iSub,rwd,r}, respPeaksDiff{iSub,rwd,r}, 1:length(resp{iSub,rwd,r}), respInterpMethod, 'extrap');
            %             interpRespPeaksDiff{iSub,rwd,r}(end:ecgRunLength) = NaN;
            %
            %             rvt{iSub,rwd,r} = (interpPeaksAmp{iSub,rwd,r} - interpTroughsAmp{iSub,rwd,r})./interpRespPeaksDiff{iSub,rwd,r};
            %
            %we will junk the first 10 TRs of BOLD anyway, so starting from
            %6 sec. BUT, this measure is lagging.
            for t=respStdWindow:length(resp{iSub,rwd,r})
                rv{iSub,rwd,r}(t) = nanstd(resp{iSub,rwd,r}(t-respStdWindow+1:t));
            end
            rv{iSub,rwd,r}(1:respStdWindow) = rv{iSub,rwd,r}(t);
            rv{iSub,rwd,r}(end:ecgRunLength) = NaN;
%             rvt{iSub,rwd,r}(end:ecgRunLength) = NaN;
            temp = reshape(rv{iSub,rwd,r}, ecgTrial,[]);
            %             rv{iSub,rwd,r}(1:respStdWindow) = squeeze(mean(temp(1:respStdWindow,:),2));%extrapolate to first trial from mean across trials
            trialsRV{iSub,rwd}(r,:,:) = reshape(rv{iSub,rwd,r}, ecgTrial,[]);
            %             trialsRVT{iSub,rwd}(r,:,:) = reshape(rvt{iSub,rwd,r}, ecgTrial,[]);
            meanTrialRvRun{iSub,rwd}(r,:) = squeeze(mean(trialsRV{iSub,rwd}(r,:,2:end),3));%mean across trials, exclude first trial
            %             meanTrialRvtRun{iSub,rwd}(r,:) = squeeze(mean(trialsRVT{iSub,rwd}(r,:,2:end),3));%mean across trials, exclude first trial
            
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
%         rwdMeanRVT(iSub,rwd,:) = mean(meanTrialRvtRun{iSub,rwd});
        rwdMeanRV(iSub,rwd,:) = mean(meanTrialRvRun{iSub,rwd});
        rwdRvTC{iSub,rwd} = [];
%         rwdRvtTC{iSub,rwd} = [];
        for r=1:numRuns(iSub,rwd)
            %removing junked frames & concatenate
            rwdRvTC{iSub,rwd} = horzcat(rwdRvTC{iSub,rwd}, rv{iSub,rwd,r}(ecgSampleRate*1.5*junkedFrames+1:end)');
%             rwdRvtTC{iSub,rwd} = horzcat(rwdRvtTC{iSub,rwd}, rvt{iSub,rwd,r}(ecgSampleRate*1.5*junkedFrames+1:end));
        end
        %downsample pulse using median
        trRV = reshape(rwdRvTC{iSub,rwd}, ecgSampleRate*1.5,[]);
%         trRVT = reshape(rwdRvtTC{iSub,rwd}, ecgSampleRate*1.5,[]);
        downsampledRV{iSub,rwd} = nanmedian(trRV);
%         downsampledRVT{iSub,rwd} = nanmedian(trRVT);
        %create design matrix for deconvolution
        designMatRv{iSub,rwd}(1,:) = downsampledRV{iSub,rwd};
%         designMatRvt{iSub,rwd}(1,:) = downsampledRVT{iSub,rwd};
        for i=2:deconvLength
            designMatRv{iSub,rwd}(i,:) = circshift(designMatRv{iSub,rwd}(i-1,:),1);
%             designMatRvt{iSub,rwd}(i,:) = circshift(designMatRvt{iSub,rwd}(i-1,:),1);
        end
%         designMatResp{iSub,rwd} = [designMatRv{iSub,rwd}; designMatRvt{iSub,rwd}];
        designMatResp{iSub,rwd} = [designMatRv{iSub,rwd}];%ONLY RV
        
        %ECG
%         rwdMeanPeaksDiff(iSub,rwd,:) = mean(meanTrialPeaksDiffRun{iSub,rwd});
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
        
        %KEEP ONLY GOOD TRIALS
        designMatPulse{iSub,rwd} = designMatPulse{iSub,rwd}(:,goodTRs);
        designMatResp{iSub,rwd} = designMatResp{iSub,rwd}(:,goodTRs);
        designMatRespPulse{iSub,rwd} = designMatRespPulse{iSub,rwd}(:,goodTRs);
        
        %compute Fourier amp and angle of task-related response, for ALL voxels
        fullFft = fft(allTseries,[],2);
        fullFft = fullFft(:,1:(1+size(fullFft,2)/2));
        fundFreq = 1+size(allTseries,2)/trialLength;
        f = fullFft(:,fundFreq);
        allVoxTaskPhase{iSub,rwd} = angle(f);
        allVoxTaskAmp{iSub,rwd} = abs(f);
        absFullFft = sum(abs(fullFft),2);
        allVoxTaskCo{iSub,rwd} = abs(f)./absFullFft;
        
        temp = reshape(allTseries, size(allTseries,1), trialLength, size(allTseries,2)/trialLength);%(vox,trialLength,numTrials)
        allVoxTrialResponse{iSub,rwd} = mean(temp,3);%mean over trials.(vox,trialLength)
        temp = fft(allVoxTrialResponse{iSub,rwd},[],2);
        f = temp(:,2);
        allVoxTaskPhase{iSub,rwd} = angle(f);
        allVoxTaskAmp{iSub,rwd} = abs(f);
        
        
        
            
        
        for iRoi = 1:length(roiNames)

            roiTC{iSub,iRoi,rwd} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
            roiTC{iSub,iRoi,rwd}.tSeries = 100*(roiTC{iSub,iRoi,rwd}.tSeries-1);%percent signal change
            if toZscore
                roiTC{iSub,iRoi,rwd}.tSeries = zscoreConcat(roiTC{iSub,iRoi,rwd}.tSeries, concatInfo{iSub,rwd});
            end
            
            if regressGlobalMean
                regressBetasGlobal{iSub,rwd,iRoi} = globalMean{iSub,rwd}\roiTC{iSub,iRoi,rwd}.tSeries';
                roiTC{iSub,iRoi,rwd}.tSeries = roiTC{iSub,iRoi,rwd}.tSeries - regressBetasGlobal{iSub,rwd,iRoi}'*globalMean{iSub,rwd}';
            end

            roiMeanTseries{iSub,iRoi,rwd}(:) = nanmean(roiTC{iSub,iRoi,rwd}.tSeries);%mean across voxels
            %             foo = nanmean(rois{iRoi}{rwd}.tSeries(goodVox{iRoi}{selectionRun}==1,:))-1;
            subTrialResponse{iSub,iRoi,rwd} = reshape(roiMeanTseries{iSub,iRoi,rwd}(:), trialLength, length(roiMeanTseries{iSub,iRoi,rwd}(:))/trialLength);
            

            
            %average per run - ALL TRIALS AVERAGED, NOT ONLY CORRECT TRIALS
            reshapedTrials = reshape(subTrialResponse{iSub,iRoi,rwd},10*15,[]);%divide into runs
            subRoiRuns{iSub,iRoi,rwd} = reshapedTrials;
            runFFT = abs(fft(reshapedTrials));
            runMeanFFT(iSub,iRoi,rwd,:) = mean(runFFT,2);
            
            subMeanRunTC(iSub,iRoi,rwd,:) = squeeze(mean(reshapedTrials,2));%average over runs
            subStdRunTC(iSub,iRoi,rwd,:) = squeeze(std(reshapedTrials,0,2));%std over runs
            
            
            subTrialResponse{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:,goodTrials==1);
            reshapedTrials = reshape(subTrialResponse{iSub,iRoi,rwd},trialLength,[]);
            if iSub==1
                allTrials{iRoi,rwd} = reshapedTrials;
            else
                allTrials{iRoi,rwd} = [allTrials{iRoi,rwd} reshapedTrials];
            end
            %average per subject
            subResponse(iSub,iRoi,rwd,:) = mean(subTrialResponse{iSub,iRoi,rwd},2);
            subStd(iSub,iRoi,rwd,:) = std(subTrialResponse{iSub,iRoi,rwd},0,2);
            subplot(rows,cols,iSub);
            plot(squeeze(subResponse(iSub,iRoi,rwd,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
            hold on
            
            
            %single voxel responses
            voxTrials{iSub,iRoi,rwd} = reshape(roiTC{iSub,iRoi,rwd}.tSeries, size(roiTC{iSub,iRoi,rwd}.tSeries,1), trialLength, length(roiMeanTseries{iSub,iRoi,rwd}(:))/trialLength);
            voxGoodTrials{iSub,iRoi,rwd} = voxTrials{iSub,iRoi,rwd}(:,:,goodTrials);
            meanVoxTrial{iSub,iRoi,rwd} = squeeze(nanmean(voxGoodTrials{iSub,iRoi,rwd},3));
            
            
        end
    end
    
    
    
    
    %load benson eccentricity maps
    v = viewSet(v, 'curGroup', 'templates');
    templateGroup = viewGet(v,'curGroup');
    v = loadAnalysis(v, 'mrDispOverlayAnal/templateRet.mat');
    for iRoi = 1:length(roiNames)
        bensonData = loadROIbensonMatching(v,roiNames{iRoi},1,templateGroup,1,concatGroupNum);
        eccen{iSub,iRoi} = bensonData{1}.eccen;
        ang{iSub,iRoi} = bensonData{1}.ang;
        areas{iSub,iRoi} = bensonData{1}.areas;
    end
    
    title(getLastDir(pwd));
    
    deleteView(v);
    cd ..
end
set(gcf,'position',[100 100 1000 500]);
%%
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
%%
save([dataFolder 'rwdTC_physio' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
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
    'rwdRvTC','downsampledRV',...
    'designMatPulse','designMatRespPulse','designMatResp','deconvLength',...
    'allGoodTrials');

runRwd
toc

%% plot polar histogram of task related phase
figure
rows=2;
cols=numSubs;
nbins=20;
nVox=1000;
thresh=0.015;
for iSub=1:numSubs
    %     bestVox = sort(allVoxTaskCo{iSub,1}+allVoxTaskCo{iSub,2},'descend');
    for rwd=1:2
        subplot(rows,cols,iSub+(rwd-1)*cols)
        %         polarhistogram(allVoxTaskPhase{iSub,rwd}(bestVox(1:nVox)),nbins);
        polarhistogram(allVoxTaskPhase{iSub,rwd}(allVoxTaskCo{iSub,rwd}>thresh),nbins);
        set(gca,'RTickLabel',[]);
        set(gca,'ThetaTickLabel',[]);
    end
end
%%
function newConcat = zscoreConcat(concatData, concatInfo)%voxels,TRs
newConcat = [];
for iRun=1:size(concatInfo.runTransition,1)
    thisRun = concatData(:,concatInfo.runTransition(iRun,1):concatInfo.runTransition(iRun,2));
    newConcat = [newConcat zscore(thisRun,0,2)];
end
% newConcat = concatData;
end