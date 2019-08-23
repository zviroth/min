onlyCorrect=0;%1=correct,2=incorrect,0=all
toZscore=1;%0 or 1
curFolder = pwd;
% dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', ...
    '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', ...
    '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621',...
    'extra_000820171219', 'extra_001620171213', '003520180212', '004020180221',...
    '004120180223'};
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', ...
%     '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', ...
%     '0042i20180412', '0045i20180309', '0046i20180409',  '005220180621',...
%     'extra_000820171219', 'extra_001620171213', '003520180212', '004020180221',...
%     '004120180223'};

numSubs=length(subFolders);
% roiNames = {'rV1_eccen8', 'lV1_eccen8','rV2_eccen8','lV2_eccen8','rV3_eccen8','lV3_eccen8','rh_17Networks_16','lh_17Networks_16'};
% roiNames = {'leftBenson', 'rightBenson'};
roiNames = { 'leftBenson', 'rightBenson','lV1_eccen8','rV1_eccen8','lV2_eccen8','rV2_eccen8','lV3_eccen8','rV3_eccen8','lh_17Networks_16','rh_17Networks_16'};


cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end
trialLength=10;
clear motionEst d concatInfo sumResponse subResponse roiMeanTseries meanResponse roiTC subRunResponse trialCorrectness trialResponse trialRT propCorrect subMeanRunTC
plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2],[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.'};
rows=2;
cols=ceil(numSubs/rows);
for iSub = 1:numSubs
    motionDir = ['/Volumes/MH02086153MACDT-Drobo/allMinSubjects/dfiles/' subFolders{iSub} '/'];
    cd(subFolders{iSub});
    v=newView;
    
%     % Find the correct Concatenation group
%     concatGroupNum = viewGet(v,'groupNum','concat'); %concatenation of multiple sessions
%     if isempty(concatGroupNum)%single session
%         concatGroupNum = viewGet(v,'groupNum','Concatenation');
%     end
    concatGroupNum = viewGet(v,'groupNum','Concatenation');
    
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', concatGroupNum);
    nScans = viewGet(v, 'nscans');
    clear rois
    for iRoi = 1:length(roiNames)
        for iScan = 1:nScans%2 concatenations, 1 for each reward type
            s = viewGet(v, 'stimfile', iScan);
            if ~isfield(s{1}, 'stimulus')
                s=s{1};
            end
%                 rwdType = s{1}.stimulus.rewardVal;
%             else
%                 rwdType = s{1}{1}.stimulus.rewardVal;
%             end
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
            
            for r=1:length(s)
%                 temp = load([motionDir 'dfile.r' concatInfo{iSub,rwd}.filename{r}(2:3) '.1D']);
%                 motionEst{iSub,rwd}(r,:,:) = temp(1:156,:);%for subject 49 some runs were longer???
               motionEst{iSub,rwd,r} = load([motionDir 'dfile.r' concatInfo{iSub,rwd}.filename{r}(2:3) '.1D']);
                
            end
            
            roiTC{iSub,iRoi,rwd} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
            roiTC{iSub,iRoi,rwd}.tSeries = 100*(roiTC{iSub,iRoi,rwd}.tSeries-1);%percent signal change
            if toZscore
                roiTC{iSub,iRoi,rwd}.tSeries = zscoreConcat(roiTC{iSub,iRoi,rwd}.tSeries, concatInfo{iSub,rwd});
            end

            for r=1:length(s)
                %if subject stopped responding at a certain point we
                %need to fill in the last incorrect trials:
                runRwd{iSub,rwd}(r) = s{r}.stimulus.rewardVal;
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
            
            
            roiMeanTseries{iSub,iRoi,rwd}(:) = nanmean(roiTC{iSub,iRoi,rwd}.tSeries);%mean across voxels
            %             foo = nanmean(rois{iRoi}{rwd}.tSeries(goodVox{iRoi}{selectionRun}==1,:))-1;
            subTrialResponse{iSub,iRoi,rwd} = reshape(roiMeanTseries{iSub,iRoi,rwd}(:), trialLength, length(roiMeanTseries{iSub,iRoi,rwd}(:))/trialLength);
            temp = trialCorrectness{iSub,rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);
            temp = trialResponse{iSub,rwd}(:,2:end-1);
            trialResponseVec = temp(:);
            

            %average per run - ALL TRIALS AVERAGED, NOT ONLY CORRECT TRIALS
            reshapedTrials = reshape(subTrialResponse{iSub,iRoi,rwd},10*15,[]);%divide into runs
            subRoiRuns{iSub,iRoi,rwd} = reshapedTrials;
            subMeanRunTC(iSub,iRoi,rwd,:) = squeeze(mean(reshapedTrials,2));%average over runs
            subStdRunTC(iSub,iRoi,rwd,:) = squeeze(std(reshapedTrials,0,2));%std over runs
            
            if onlyCorrect ==1 %ONLY CORRECT
                goodTrials = trialCorrectnessVec==1;
                %                 numTrials = sum(trialCorrectnessVec);
            elseif onlyCorrect ==2 % ONLY INCORRECT!!!
                goodTrials = trialCorrectnessVec==0 & trialResponseVec>0;
            else % including all trials with a response
                goodTrials = trialResponseVec>0;
            end
            subTrialResponse{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:,goodTrials);
 
            
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
end
zScoreString = '';
if toZscore
    zScoreString = '_zscored';
end
%%
save([dataFolder 'rwdTC_motion' onlyCorrectString zScoreString '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'motionEst');



%%
function newConcat = zscoreConcat(concatData, concatInfo)%voxels,TRs
newConcat = [];
for iRun=1:size(concatInfo.runTransition,1)
    thisRun = concatData(:,concatInfo.runTransition(iRun,1):concatInfo.runTransition(iRun,2));
    newConcat = [newConcat zscore(thisRun,0,2)];
end
% newConcat = concatData;
end