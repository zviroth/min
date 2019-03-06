
roiNames = {'rV1_eccen8', 'rh_17Networks_16'};
cd '/Volumes/MH02086153MACDT-Drobo/decodingAnalysis/rwd/';
subdirs = dir('s00*');
% subdirs = subdirs(1:3);
% subdirs = dir('s004020180221');
cols=ceil(length(subdirs)/2);
rows = ceil(length(subdirs)/cols);
mrQuit;
if exist('v')
    deleteView(v);
end
plotColors = {[0 0 1], [1 0 0]};
plotStyles = {'-','--'};
clear d concatInfo sumResponse subResponse roiMeanTseries meanResponse rois
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
                rois{iRoi,1} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
                concatInfo{iSub,1} = viewGet(v, 'concatInfo', iScan);
            elseif strcmp(rwdType, 'L')
                rois{iRoi,2} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
                concatInfo{iSub,2} = viewGet(v, 'concatInfo', iScan);
            else
                disp('wtf');
                keyboard
            end
        end
        for rwd=1:2

            roiMeanTseries{iSub,iRoi,rwd}(:) = nanmean(rois{iRoi,rwd}.tSeries);
            %             foo = nanmean(rois{iRoi}{rwd}.tSeries(goodVox{iRoi}{selectionRun}==1,:))-1;
            subResponse{iSub,iRoi}(rwd,:) = mean(reshape(roiMeanTseries{iSub,iRoi,rwd}(:), 10, length(roiMeanTseries{iSub,iRoi,rwd}(:))/10),2);
            if iSub==1
                sumResponse{rwd,iRoi} = subResponse{iSub,iRoi}(rwd,:);
            else
                sumResponse{rwd,iRoi} = sumResponse{rwd,iRoi}+subResponse{iSub,iRoi}(rwd,:);
            end
            plot(subResponse{iSub,iRoi}(rwd,:), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
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
        
        meanResponse{rwd,iRoi} = sumResponse{rwd,iRoi}./length(subdirs);
        plot(meanResponse{rwd,iRoi}(:), plotStyles{iRoi}, 'Color', plotColors{rwd}, 'linewidth', 1);
        hold on
    end
end
legend('high V1', 'high DMN', 'low V1', 'low DMN');

% %%
% figure(3)
% plot(corrThisSub)

%% make response templates using all runs and both reward levels, and measure response amplitude in each run and each reward level separately
ntrials=15;
for iSub = 1:length(subdirs)
    for iRoi= 1:length(roiNames)
        subMeanResponse{iSub,iRoi} = mean(subResponse{iSub,iRoi}(:,:));%averaged over reward type
        for rwd=1:2
            for iRun = 1:concatInfo{iSub,rwd}.n
                runTseries = roiMeanTseries{iSub,iRoi,rwd}(concatInfo{iSub,rwd}.runTransition(iRun,1):concatInfo{iSub,rwd}.runTransition(iRun,2));
%we can either get an amplitude for each trial or for each run
                runMean = mean(reshape(runTseries, 10, length(runTseries)/10),2);
                amp = runMean\subMeanResponse{iSub,iRoi}';
%                 [c, lag, amp] = myxcorr(runTseries,subMeanResponse{iSub,iRoi});
                
                runAmp{iSub,rwd}(iRoi,iRun) = amp;
                runLag{iSub,rwd}(iRoi,iRun) = lag;
                
            end
            meanAmp(iSub,rwd,iRoi) = mean(runAmp{iSub,rwd}(iRoi,:));
%             meanLag(iSub,rwd,iRoi) = mean(runLag{iSub,rwd}(iRoi,:));
        end
    end
end