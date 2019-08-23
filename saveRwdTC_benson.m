onlyCorrect=0;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
numSubs=length(subFolders);
roiNames = {'rV1_eccen8', 'lV1_eccen8','rV2_eccen8','lV2_eccen8','rV3_eccen8','lV3_eccen8','rh_17Networks_16','lh_17Networks_16'};
roiNames = {'leftBenson', 'rightBenson'};

cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end
trialLength=10;

for iSub = 1:numSubs
    if strfind(subdirs(iSub).name, 'i')
        hasEyeData(iSub) = 1;
    else
        hasEyeData(iSub) = 0;
    end
    cd(subFolders{iSub});
    v=newView;
    
    % Find the correct Concatenation group
    concatGroupNum = viewGet(v,'groupNum','concat'); %concatenation of multiple sessions
    if isempty(concatGroupNum)%single session
        concatGroupNum = viewGet(v,'groupNum','Concatenation');
    end
    
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', 'Concatenation');
    nScans = viewGet(v, 'nscans');
    
end


% subdirs=[];
% for i=1:length(dataFolder)
%     cd(dataFolder{i});
%     temp = dir('s00*');
%     if length(temp)==0
%         temp = dir('00*');
%     end
%     subdirs = [subdirs; temp];
% end

% subdirs = subdirs(1:3);
% subdirs = dir('s004020180221');
cols=ceil(length(subdirs)/2);
rows = ceil(length(subdirs)/cols);
plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
clear d concatInfo sumResponse subResponse roiMeanTseries meanResponse roiTC subRunResponse trialCorrectness trialResponse trialRT propCorrect


for iSub = 1:length(subdirs)
    subplot(rows,cols,iSub);
    % go into the subject's directory
    disp(sprintf('Processing %s ...', subdirs(iSub).name));
    cd(fullfile(subdirs(iSub).folder, subdirs(iSub).name))
    
    v=newView;
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', 'Concatenation');
    concatGroupNum = viewGet(v,'curGroup');
    % load the event related analysis analysis
%     v = loadAnalysis(v, 'glmAnalStats/GLM.mat');
    
    % load the data
    nScans = viewGet(v, 'nscans');
    
    clear rois
    for iRoi = 1:length(roiNames)
        for iScan = 1:nScans%2 concatenations, 1 for each reward type
%             d{iScan} = viewGet(v, 'd', iScan);
            %             concatInfo{iSub,iScan} = viewGet(v, 'concatInfo', iScan);
            s = viewGet(v, 'stimfile', iScan);
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
            roiTC{iSub,iRoi,rwd} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
            roiTC{iSub,iRoi,rwd}.tSeries = zscoreConcat(roiTC{iSub,iRoi,rwd}.tSeries, concatInfo{iSub,rwd});
            
            
            
            for r=1:length(s)
                    %if subject stopped responding at a certain point we
                    %need to fill in the last incorrect trials:
                    trialCorrectness{iSub,rwd}(r,:) = [s{r}.stimulus.trialCorrectness zeros(1,17-size(s{r}.stimulus.trialCorrectness,2))];
%                     trialResponse{iSub,rwdTypeNum,r} = s{r}.stimulus.trialResponse;
                    trialResponse{iSub,rwd}(r,:) = [s{r}.stimulus.trialResponse zeros(1,17-size(s{r}.stimulus.trialResponse,2))];%0 means no response

                    trialRT{iSub,rwd,r} = s{r}.stimulus.trialRT;
                    propCorrect{iSub,rwd}(r) = s{r}.stimulus.percentCorrect;
                    stairThresh{iSub,rwd}(r,:) = [s{r}.stimulus.stair{1}.threshold s{r}.stimulus.stair{2}.threshold];
                    
            end
            expName{iSub,rwd} = s{iScan}.task{1}.taskFilename;
            
            %if there is eye data, extract it
            if hasEyeData(iSub)
                for r=1:length(s)
                    stimfile = s{r}.filename;
                    e{r} = myGetTaskEyeTraces(stimfile, 'removeBlink=3');
                    runEyeSize{rwd}(r,:) = size(e{r}.eye.pupil);
                end
                
                trialLengthEye(rwd) = max(runEyeSize{rwd}(:,2));
                numTrials(rwd) = sum(runEyeSize{rwd}(:,1));
                rwdPupil{iSub,rwd} = NaN(numTrials(rwd), trialLengthEye(rwd));
                trialCounter=0;
                for r=1:length(s)
                    rwdPupil{iSub,rwd}(trialCounter+1:trialCounter+runEyeSize{rwd}(r,1), 1:runEyeSize{rwd}(r,2)) = e{r}.eye.pupil;
                    trialCounter = trialCounter + runEyeSize{rwd}(r,1);
                end
                meanPupil{iSub,rwd} = nanmean(rwdPupil{iSub,rwd})';
            else
                meanPupil{iSub,rwd} = [];
                rwdPupil{iSub,rwd}=[];
            end
            
        end
        for rwd=1:2
            roiMeanTseries{iSub,iRoi,rwd}(:) = nanmean(roiTC{iSub,iRoi,rwd}.tSeries);%mean across voxels
            %             foo = nanmean(rois{iRoi}{rwd}.tSeries(goodVox{iRoi}{selectionRun}==1,:))-1;
            subTrialResponse{iSub,iRoi,rwd} = reshape(roiMeanTseries{iSub,iRoi,rwd}(:), trialLength, length(roiMeanTseries{iSub,iRoi,rwd}(:))/trialLength);
            temp = trialCorrectness{iSub,rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);
            temp = trialResponse{iSub,rwd}(:,2:end-1);
            trialResponseVec = temp(:);
            
            %             numTrials=length(trialCorrectnessVec);
            
            %average per run - ALL TRIALS AVERAGED, NOT ONLY CORRECT TRIALS
            reshapedTrials = reshape(subTrialResponse{iSub,iRoi,rwd},10,15,[]);
            subRunResponse{iSub,iRoi,rwd} = squeeze(mean(reshapedTrials,2));
            
            
            if onlyCorrect ==1 %ONLY CORRECT
                goodTrials = trialCorrectnessVec==1;
                %                 numTrials = sum(trialCorrectnessVec);
            elseif onlyCorrect ==2 % ONLY INCORRECT!!!
                goodTrials = trialCorrectnessVec==0 & trialResponseVec>0;
            else % including all trials with a response
                goodTrials = trialResponseVec>0;
            end
            subTrialResponse{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:,goodTrials);
            
%             if onlyCorrect ==1 %ONLY CORRECT
%                 subTrialResponse{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:,trialCorrectnessVec==1);
%                 %                 numTrials = sum(trialCorrectnessVec);
%             elseif onlyCorrect ==2 % ONLY INCORRECT!!!
%                 subTrialResponse{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:,trialCorrectnessVec==0 & trialResponseVec>0);
%             else % including all trials with a response
%                 subTrialResponse{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:,trialResponseVec>0);
%             end
            
            reshapedTrials = reshape(subTrialResponse{iSub,iRoi,rwd},trialLength,[]);
            if iSub==1
                allTrials{iRoi,rwd} = reshapedTrials;
            else
                allTrials{iRoi,rwd} = [allTrials{iRoi,rwd} reshapedTrials];
            end
            %average per subject
            subResponse(iSub,iRoi,rwd,:) = mean(subTrialResponse{iSub,iRoi,rwd},2);

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
end
set(gcf,'position',[100 100 1000 500]);
%%
figure(2)
for rwd=1:2
    for iRoi = 1:length(roiNames)
        
%         meanResponse(rwd,iRoi,:) = sumResponse(rwd,iRoi,:)./length(subdirs);
        meanResponse(rwd,iRoi,:) = mean(subResponse(:,iRoi,rwd,:));
        stdResponse(rwd,iRoi,:) = std(subResponse(:,iRoi,rwd,:));
        plot(squeeze(meanResponse(rwd,iRoi,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
        hold on
    end
end

onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
end
save([saveFolder 'rwdTC_benson' onlyCorrectString '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', 'meanResponse', 'stdResponse',...
    'roiTC', 'allTrials', ...
    'subdirs', 'roiNames','subTrialResponse','subRunResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'rwdPupil','meanPupil','hasEyeData','expName','stairThresh','eccen','ang','areas','trialLength');



%%
function newConcat = zscoreConcat(concatData, concatInfo)%voxels,TRs
% newConcat = [];
% for iRun=1:size(concatInfo.runTransition,1)
%     thisRun = concatData(:,concatInfo.runTransition(iRun,1):concatInfo.runTransition(iRun,2));
%     newConcat = [newConcat zscore(thisRun,0,2)];
% end
newConcat = concatData;
end