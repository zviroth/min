tic

onlyCorrect=0;
clear dataFolder
dataFolder{1} = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects/';

roiNames = {'rV1_eccen8', 'lV1_eccen8','rV2_eccen8','lV2_eccen8','rV3_eccen8','lV3_eccen8','rh_17Networks_16','lh_17Networks_16'};
roiNames = {'rV1_eccen8', 'lV1_eccen8'};
subdirs=[];
for i=1:length(dataFolder)
    cd(dataFolder{i});
    temp = dir('s00*');
    if length(temp)==0
        temp = dir('00*');
    end
    subdirs = [subdirs; temp];
end

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
    subplot(rows,cols,iSub);
    % go into the subject's directory
    disp(sprintf('Processing %s ...', subdirs(iSub).name));
    cd(fullfile(subdirs(iSub).folder, subdirs(iSub).name))
    
    v=newView;
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', 'Concatenation');
    
%     % load the event related analysis analysis
%     v = loadAnalysis(v, 'glmAnalStats/GLM.mat');
%     
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
            
            roiTC{iRoi,rwdTypeNum} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);

            for r=1:length(s)
                    %if subject stopped responding at a certain point we
                    %need to fill in the last incorrect trials:
                    trialCorrectness{iSub,rwdTypeNum}(r,:) = [s{r}.stimulus.trialCorrectness zeros(1,17-size(s{r}.stimulus.trialCorrectness,2))];
            end

            
        end
        for rwd=1:2
            roiMeanTseries{iSub,iRoi,rwd}(:) = nanmean(roiTC{iRoi,rwd}.tSeries);%mean across voxels
            %             foo = nanmean(rois{iRoi}{rwd}.tSeries(goodVox{iRoi}{selectionRun}==1,:))-1;
            subTrialResponse{iSub,iRoi,rwd} = reshape(roiMeanTseries{iSub,iRoi,rwd}(:), 10, length(roiMeanTseries{iSub,iRoi,rwd}(:))/10);
            temp = trialCorrectness{iSub,rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);

            if onlyCorrect
                subTrialResponse{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:,trialCorrectnessVec==1);
                %                 numTrials = sum(trialCorrectnessVec);
            end

            %average per subject
            subResponse(iSub,iRoi,rwd,:) = mean(subTrialResponse{iSub,iRoi,rwd},2);

            plot(1:1.5:15,100*squeeze(subResponse(iSub,iRoi,rwd,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
            hold on
        end
    end
    
    title(getLastDir(pwd));
    
    deleteView(v);
end
set(gcf,'position',[100 100 1000 500]);
xlabel('time (s)');
ylabel('BOLD response (% signal change)');
%%
figure(2)
clf
rows=1;
cols=length(roiNames);
for iRoi = 1:length(roiNames)
    subplot(rows,cols,iRoi)
    for rwd=1:2
%         meanResponse(rwd,iRoi,:) = sumResponse(rwd,iRoi,:)./length(subdirs);
        meanResponse(rwd,iRoi,:) = mean(subResponse(:,iRoi,rwd,:));
        stdResponse(rwd,iRoi,:) = std(subResponse(:,iRoi,rwd,:));
        plot(squeeze(meanResponse(rwd,iRoi,:)), 'Color', plotColors{rwd}, 'linewidth', 1);
        hold on
    end
    title(roiNames{iRoi})
end

toc