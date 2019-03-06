onlyCorrect=0;
roiNames = {'rV1_eccen8', 'rh_17Networks_16','lV1_eccen8.mat','lh_17Networks_16.mat'};
roiNames = {'rV1_eccen8', 'lV1_eccen8.mat','rV1_localizer','rV2_eccen8','lV2_eccen8','rh_17Networks_16','lh_17Networks_16.mat'};
roiNames = {'Rviz_localizer','Lviz_localizer','rV1_eccen8', 'lV1_eccen8.mat','rV2_eccen8','lV2_eccen8','rh_17Networks_16','lh_17Networks_16.mat'};
% roiNames = {'rV1_eccen8'};
dataFolder = '/Volumes/MH02086153MACDT-Drobo/decodingAnalysis/rwd/';
% eval(['cd ' dataFolder]);
cd(dataFolder);
subdirs = dir('s00*');
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
    subplot(rows,cols,iSub);
    % go into the subject's directory
    disp(sprintf('Processing %s ...', subdirs(iSub).name));
    cd(fullfile(subdirs(iSub).folder, subdirs(iSub).name))
    
    v=newView;
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', 'Concatenation');
    
    % load the event related analysis analysis
    v = loadAnalysis(v, 'glmAnalStats/GLM.mat');
    
    % load the data
    nScans = viewGet(v, 'nscans');
    
    clear rois
    for iRoi = 1:length(roiNames)
        for iScan = 1:nScans%2 concatenations, 1 for each reward type
            d{iScan} = viewGet(v, 'd', iScan);
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
            concatInfo{iSub,rwdTypeNum} = viewGet(v, 'concatInfo', iScan);
            
            for r=1:length(s)
                    %if subject stopped responding at a certain point we
                    %need to fill in the last incorrect trials:
                    trialCorrectness{iSub,rwdTypeNum}(r,:) = [s{r}.stimulus.trialCorrectness zeros(1,17-size(s{r}.stimulus.trialCorrectness,2))];
                    trialResponse{iSub,rwdTypeNum,r} = s{r}.stimulus.trialResponse;
                    trialRT{iSub,rwdTypeNum,r} = s{r}.stimulus.trialRT;
                    propCorrect{iSub,rwdTypeNum}(r) = s{r}.stimulus.percentCorrect;

            end
            
        end
        for rwd=1:2
            roiMeanTseries{iSub,iRoi,rwd}(:) = nanmean(roiTC{iRoi,rwd}.tSeries);%mean across voxels
            %             foo = nanmean(rois{iRoi}{rwd}.tSeries(goodVox{iRoi}{selectionRun}==1,:))-1;
            subTrialResponse{iSub,iRoi,rwd} = reshape(roiMeanTseries{iSub,iRoi,rwd}(:), 10, length(roiMeanTseries{iSub,iRoi,rwd}(:))/10);
            temp = trialCorrectness{iSub,rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);
            %             numTrials=length(trialCorrectnessVec);
            
            %average per run - ALL TRIALS AVERAGED, NOT ONLY CORRECT TRIALS
            reshapedTrials = reshape(subTrialResponse{iSub,iRoi,rwd},10,15,[]);
            subRunResponse{iSub,iRoi,rwd} = squeeze(mean(reshapedTrials,2));
            
            if onlyCorrect
                subTrialResponse{iSub,iRoi,rwd} = subTrialResponse{iSub,iRoi,rwd}(:,trialCorrectnessVec==1);
                %                 numTrials = sum(trialCorrectnessVec);
            end
            reshapedTrials = reshape(subTrialResponse{iSub,iRoi,rwd},10,[]);
            if iSub==1
                allTrials{iRoi,rwd} = reshapedTrials;
            else
                allTrials{iRoi,rwd} = [allTrials{iRoi,rwd} reshapedTrials];
            end
            %average per subject
            subResponse(iSub,iRoi,rwd,:) = mean(subTrialResponse{iSub,iRoi,rwd},2);
%             if iSub==1
%                 sumResponse(rwd,iRoi,:) = squeeze(subResponse(iSub,iRoi,rwd,:));
%             else
%                 sumResponse(rwd,iRoi,:) = sumResponse(rwd,iRoi,:)+squeeze(subResponse(iSub,iRoi,rwd,:));
%             end
            plot(squeeze(subResponse(iSub,iRoi,rwd,:)), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
            hold on
        end
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




% legend('high right V1', 'high right DMN', 'high left V1', 'high left DMN', 'low right V1', 'low right DMN', 'low left V1','low left DMN');

% legend('high right V1', 'low right V1','high right DMN', 'low right DMN', 'high left V1', 'low left V1', 'high left DMN','low left DMN');

% %%
% figure(3)
% plot(corrThisSub)
onlyCorrectString = '';
if onlyCorrect
    onlyCorrectString = '_correct';
end
save([dataFolder 'subMeanRoiTC' onlyCorrectString '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', 'meanResponse', 'stdResponse',...
    'roiTC', 'allTrials', ...
    'subdirs', 'roiNames','subTrialResponse','subRunResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect');

% %% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
% ntrials=15;
% for iSub = 1:length(subdirs)
%     for iRoi= 1:length(roiNames)
%         subMeanResponse{iSub,iRoi} = mean(subResponse{iSub,iRoi}(:,:));%averaged over reward type
%         for rwd=1:2
%             for iRun = 1:concatInfo{iSub,rwd}.n
%                 runTseries = roiMeanTseries{iSub,iRoi,rwd}(concatInfo{iSub,rwd}.runTransition(iRun,1):concatInfo{iSub,rwd}.runTransition(iRun,2));
% %we can either get an amplitude for each trial or for each run
%                 runMean = mean(reshape(runTseries, 10, length(runTseries)/10),2);
%                 amp = runMean\subMeanResponse{iSub,iRoi}';
% %                 [c, lag, amp] = myxcorr(runTseries,subMeanResponse{iSub,iRoi});
%
%                 runAmp{iSub,rwd}(iRoi,iRun) = amp;
%                 runLag{iSub,rwd}(iRoi,iRun) = lag;
%
%             end
%             meanAmp(iSub,rwd,iRoi) = mean(runAmp{iSub,rwd}(iRoi,:));
% %             meanLag(iSub,rwd,iRoi) = mean(runLag{iSub,rwd}(iRoi,:));
%         end
%     end
% end