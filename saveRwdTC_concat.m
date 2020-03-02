close all
clear all
tic
maxRT=4000;
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean = 0;%1=yes, 0=no, 2=only the Benson ROI.
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', ...
    '003220180105', '0034i20180209', '003520180328', '0039i20180220', '004020180328','004120180320', ...
    '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', ...
%     '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', ...
%     '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
numSubs=length(subFolders);
% roiNames = {'rV1_eccen8', 'lV1_eccen8','rV2_eccen8','lV2_eccen8','rV3_eccen8','lV3_eccen8','rh_17Networks_16','lh_17Networks_16'};
% roiNames = {'leftBenson', 'rightBenson'};
% roiNames = { 'leftBenson', 'rightBenson','lV1_eccen8','rV1_eccen8','lV2_eccen8','rV2_eccen8','lV3_eccen8','rV3_eccen8','lh_17Networks_16','rh_17Networks_16'};
roiNames = { 'leftBenson', 'rightBenson','Lviz_localizer','Rviz_localizer','lh_17Networks_16','rh_17Networks_16','leftCerebellarCortex','rightCerebellarCortex'};


cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end
trialLength=10;
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
        concatInfo{iSub,rwd} = viewGet(v, 'concatInfo', iScan);
        
        %get global mean
        ts = loadTSeries(v,iScan);
        allTseries = reshape(ts,[],size(ts,4));%voxels,TRs
        
%         roiTC{iSub,iRoi,rwd} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
%         roiTC{iSub,iRoi,rwd}.tSeries = 100*(roiTC{iSub,iRoi,rwd}.tSeries-1);%percent signal change
        
        benson1 = loadROITSeries(v, roiNames{1}, iScan, [], 'keepNAN',true);
        benson2 = loadROITSeries(v, roiNames{2}, iScan, [], 'keepNAN',true);
        visualRoi = 100*[benson1.tSeries-1; benson2.tSeries-1];%percent signal change
        
        if toZscore
            allTseries = zscoreConcat(allTseries, concatInfo{iSub,rwd});
            visualRoi= zscoreConcat(visualRoi, concatInfo{iSub,rwd});
        end
        globalMean{iSub,rwd} = nanmean(allTseries)';%(vox,T)
        visualMean{iSub,rwd} = nanmean(visualRoi)';%(vox,T)
        
        
        
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
        
        for iRoi = 1:length(roiNames)
            

            roiTC{iSub,iRoi,rwd} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
            roiTC{iSub,iRoi,rwd}.tSeries = 100*(roiTC{iSub,iRoi,rwd}.tSeries-1);%percent signal change
            if toZscore
               roiTC{iSub,iRoi,rwd}.tSeries = zscoreConcat(roiTC{iSub,iRoi,rwd}.tSeries, concatInfo{iSub,rwd});
            end

            if regressGlobalMean==1
                regressBetasGlobal{iSub,rwd,iRoi} = globalMean{iSub,rwd}\roiTC{iSub,iRoi,rwd}.tSeries';
                roiTC{iSub,iRoi,rwd}.tSeries = roiTC{iSub,iRoi,rwd}.tSeries - regressBetasGlobal{iSub,rwd,iRoi}'*globalMean{iSub,rwd}';
            elseif regressGlobalMean==2
                regressBetasGlobal{iSub,rwd,iRoi} = visualMean{iSub,rwd}\roiTC{iSub,iRoi,rwd}.tSeries';
                roiTC{iSub,iRoi,rwd}.tSeries = roiTC{iSub,iRoi,rwd}.tSeries - regressBetasGlobal{iSub,rwd,iRoi}'*visualMean{iSub,rwd}';
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
if regressGlobalMean==1
    globalMeanString = '_globalRegressed';
elseif regressGlobalMean==2
    globalMeanString = '_bensonRegressed';
end

ConcatProjStr = '';
if ConcatProj
    ConcatProjStr = 'ConcatProj';
end
%%
save([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'globalMean','regressBetasGlobal','runRwd',...
    'subRoiRuns','runMeanFFT',...
    'allVoxTrialResponse','allVoxTaskPhase','allVoxTaskAmp','allVoxTaskCo',...
    'voxTrials','voxGoodTrials','meanVoxTrial',...
    'maxRT','-v7.3');
%'-v7.3'
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