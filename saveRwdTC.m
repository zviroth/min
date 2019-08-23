onlyCorrect=1;
saveFolder = '/Volumes/MH02086153MACDT-Drobo/decodingAnalysis/rwd/';

% dataFolder{1} = '/Volumes/MH02086153MACDT-Drobo/decodingAnalysis/rwd/';
% roiNames = {'Rviz_localizer','Lviz_localizer','rV1_eccen8', 'lV1_eccen8.mat','rV2_eccen8','lV2_eccen8','rh_17Networks_16','lh_17Networks_16.mat'};

clear dataFolder
dataFolder{1} = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects/';
% dataFolder{1} = '/Volumes/MH02086153MACDT-Drobo/CNaP_drive_2/ryoom/data/probRwdExp/';
% dataFolder{2} = '/Volumes/MH02086153MACDT-Drobo/CNaP_drive_2/ryoom/data/AllRwds/deterministic/';
% dataFolder{3} = '/Volumes/MH02086153MACDT-Drobo/CNaP_drive_2/ryoom/data/AllRwds/probabilistic/';

roiNames = {'rV1_eccen8', 'lV1_eccen8','rV2_eccen8','lV2_eccen8','rV3_eccen8','lV3_eccen8','rh_17Networks_16','lh_17Networks_16'};
% roiNames = {'rV1_eccen8', 'lV1_eccen8.mat','rV2_eccen8','lV2_eccen8','rh_17Networks_16'};
% roiNames = {'rV1_eccen8', 'rh_17Networks_16'};

% dataFolder = '/Volumes/MH02086153MACDT-Drobo/CNaP_drive_1/RewardData/probRwdExp/';
trialLength=10;
subdirs=[];
for i=1:length(dataFolder)
    cd(dataFolder{i});
    temp = dir('s00*');
    if length(temp)==0
        temp = dir('00*');
    end
    subdirs = [subdirs; temp];
end


% subdirs = subdirs(1:3);
% subdirs = dir('s004020180221');
cols=ceil(length(subdirs)/2);
rows = ceil(length(subdirs)/cols);
mrQuit;
if exist('v')
    deleteView(v);
end
plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
clear d concatInfo sumResponse subResponse roiMeanTseries meanResponse roiTC subRunResponse trialCorrectness trialResponse trialRT propCorrect
for iSub = 1:length(subdirs)
    if strfind(subdirs(iSub).name, 'i')
        hasEyeData(iSub) = 1;
    else
        hasEyeData(iSub) = 0;
    end
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
                rwdTypeNum = 1;
            elseif strcmp(rwdType, 'L')
                rwdTypeNum = 2;
            else
                disp('wtf');
                keyboard
            end
            
            roiTC{iSub,iRoi,rwdTypeNum} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
            concatInfo{iSub,rwdTypeNum} = viewGet(v, 'concatInfo', iScan);
            
            for r=1:length(s)
                    %if subject stopped responding at a certain point we
                    %need to fill in the last incorrect trials:
                    trialCorrectness{iSub,rwdTypeNum}(r,:) = [s{r}.stimulus.trialCorrectness zeros(1,17-size(s{r}.stimulus.trialCorrectness,2))];
                    trialResponse{iSub,rwdTypeNum}(r,:) = [s{r}.stimulus.trialResponse zeros(1,17-size(s{r}.stimulus.trialResponse,2))];%0 means no response
                    trialRT{iSub,rwdTypeNum,r} = s{r}.stimulus.trialRT;
                    propCorrect{iSub,rwdTypeNum}(r) = s{r}.stimulus.percentCorrect;
                    stairThresh{iSub,rwdTypeNum}(r,:) = [s{r}.stimulus.stair{1}.threshold s{r}.stimulus.stair{2}.threshold];
                    
            end
            expName{iSub,rwdTypeNum} = s{iScan}.task{1}.taskFilename;
            
            %if there is eye data, extract it
            if hasEyeData(iSub)
                for r=1:length(s)
                    stimfile = s{r}.filename;
                    e{r} = myGetTaskEyeTraces(stimfile, 'removeBlink=3');
                    runEyeSize{rwdTypeNum}(r,:) = size(e{r}.eye.pupil);
                end
                
                trialLengthEye(rwdTypeNum) = max(runEyeSize{rwdTypeNum}(:,2));
                numTrials(rwdTypeNum) = sum(runEyeSize{rwdTypeNum}(:,1));
                rwdPupil{iSub,rwdTypeNum} = NaN(numTrials(rwdTypeNum), trialLengthEye(rwdTypeNum));
                trialCounter=0;
                for r=1:length(s)
                    rwdPupil{iSub,rwdTypeNum}(trialCounter+1:trialCounter+runEyeSize{rwdTypeNum}(r,1), 1:runEyeSize{rwdTypeNum}(r,2)) = e{r}.eye.pupil;
                    trialCounter = trialCounter + runEyeSize{rwdTypeNum}(r,1);
                end
                meanPupil{iSub,rwdTypeNum} = nanmean(rwdPupil{iSub,rwdTypeNum})';
            else
                meanPupil{iSub,rwdTypeNum} = [];
                rwdPupil{iSub,rwdTypeNum}=[];
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
            
            
            %including only correct trials
            if onlyCorrect
                subTrialResponse{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:,trialCorrectnessVec==1);
                %                 numTrials = sum(trialCorrectnessVec);
            else % including only trials with a response
                subTrialResponse{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:,trialResponseVec==1);
            end
            
            
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
if onlyCorrect
    onlyCorrectString = '_correct';
end
save([saveFolder 'rwdTC' onlyCorrectString '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', 'meanResponse', 'stdResponse',...
    'roiTC', 'allTrials', ...
    'subdirs', 'roiNames','subTrialResponse','subRunResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'rwdPupil','meanPupil','hasEyeData','expName','stairThresh','eccen','ang','areas','trialLength');

