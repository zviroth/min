tic
clear all
close all
saveFolder = '/Volumes/MH02086153MACDT-Drobo/allMinBehavioral/';
% subFolders = {'s000520180126','s000820171116','s001620171103','s002220171122','s003220180105','s003920180220',...
%     's004020180221','s004120180223','s004320180306','s004520180308','s004620180309','s004620180323','s004720180326',...
%     's004920180404'};
subFolders = {'s000520180126','s000820171116','s001620171103','s002220171122','s003220180105','s003920180220',...
    's004020180221','s004120180223','s004320180306','s004520180308','s004620180323','s004720180326',...
    's004920180404'};
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
maxRT=4000;
% for i=1:length(dataFolder)
%     cd(dataFolder{i});
%     temp = dir('s00*');
%     if length(temp)==0
%         temp = dir('00*');
%     end
%     subdirs = [subdirs; temp];
% end

for iSub = 1:length(subFolders)
    cd([saveFolder subFolders{iSub}]);
%     cd(fullfile(subdirs(iSub).folder, subdirs(iSub).name))
    stimfiles = dir('*.mat');
    for iRun=1:length(stimfiles)
        stimfile=stimfiles(iRun).name;
        s{iRun}=load(stimfile);
        
        expName{iSub,iRun} = s{iRun}.task{1}.taskFilename;
        rwdType = s{iRun}.stimulus.rewardVal;
        if strcmp(rwdType, 'H')
            rwdLevel(iSub,iRun) = 1;
        elseif strcmp(rwdType, 'L')
            rwdLevel(iSub,iRun) = 2;
        end        
        e{iRun} = myGetTaskEyeTraces(stimfile,'removeBlink=3');
        subPupil{iSub,iRun} = e{iRun}.eye.pupil;
        runSize{iSub}(iRun,:) = size(e{iRun}.eye.pupil);
    end

    
    for rwd=1:2
        rwdRuns=find(rwdLevel(iSub,:)==rwd);
        
        numRuns(rwd) = length(rwdRuns);
        numTrials(iSub,rwd) = sum(runSize{iSub}(rwdRuns,1));
        trialLength(rwd) = max(runSize{iSub}(rwdRuns,2));
        rwdPupil{iSub,rwd} = NaN(numTrials(iSub,rwd), trialLength(rwd));
        trialCounter=0;
        for irun=1:length(rwdRuns)
            r=rwdRuns(irun);
            rwdPupil{iSub,rwd}(trialCounter+1:trialCounter+runSize{iSub}(r,1), 1:runSize{iSub}(r,2)) = e{r}.eye.pupil;
            trialCounter = trialCounter + runSize{iSub}(r,1);
            
            trialCorrectness{iSub,rwd}(irun,:) = [s{r}.stimulus.trialCorrectness zeros(1,17-size(s{r}.stimulus.trialCorrectness,2))];%0 means incorrect or no response
            trialResponse{iSub,rwd}(irun,:) = [s{r}.stimulus.trialResponse zeros(1,17-size(s{r}.stimulus.trialResponse,2))];%0 means no response
            trialRT{iSub,rwd}(irun,:) = [s{r}.stimulus.trialRT zeros(1,17-size(s{r}.stimulus.trialRT,2))];%zeros means no response
            propCorrect{iSub,rwd}(irun) = s{r}.stimulus.percentCorrect;
            stairThresh{iSub,rwd}(irun,:) = [s{r}.stimulus.stair{1}.threshold s{r}.stimulus.stair{2}.threshold];
        end
        
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
        
%         if onlyCorrect ==1 %ONLY CORRECT
%                 goodTrials = trialCorrectnessVec==1;
%                 %                 numTrials = sum(trialCorrectnessVec);
%             elseif onlyCorrect ==2 % ONLY INCORRECT!!!
%                 goodTrials = trialCorrectnessVec==0 & trialResponseVec>0;
%             else % including all trials with a response
%                 goodTrials = trialResponseVec>0;
%         end
        rwdPupil{iSub,rwd} = rwdPupil{iSub,rwd}(goodTrials,:);
        subMeanCorrectness(iSub,rwd) = mean(trialCorrectness{iSub,rwd}(trialCorrectness{iSub,rwd}>0));
        trialRT{iSub,rwd} = trialRT{iSub,rwd}(goodTrials);
        subMeanRT(iSub,rwd) = mean(trialRT{iSub,rwd}(:));
        subMedianRT(iSub,rwd) = median(trialRT{iSub,rwd}(:));
        temp = mean(stairThresh{iSub,rwd});%mean across runs
        subMeanThresh(iSub,rwd) =mean(temp);%mean across staircases
        subRTvar(iSub,rwd) = std(trialRT{iSub,rwd}(:));
        
        meanPupil{iSub,rwd} = nanmean(rwdPupil{iSub,rwd})';
        stdRwd{iSub,rwd} = std(rwdPupil{iSub,rwd},0,1,'omitnan');
        stdRwd{iSub,rwd}(isnan(stdRwd{iSub,rwd})) = ones;
        meanPupil{iSub,rwd}(isnan(meanPupil{iSub,rwd}))=zeros;
        
        numTrials(iSub,rwd) = length(goodTrials);

                    
    end
    minlength = min(length(meanPupil{iSub,1}), length(meanPupil{iSub,2}));
    diffPupil{iSub} = meanPupil{iSub,1}(1:minlength) - meanPupil{iSub,2}(1:minlength);

end
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
elseif onlyCorrect==0
    onlyCorrectString = '_validresponse';
end

save([saveFolder 'behavioralData' onlyCorrectString '.mat'], 'subFolders', 'subPupil', 'runSize', 'rwdPupil',...
    'meanPupil','stdRwd','diffPupil','expName',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh','subRTvar','rwdLevel','numTrials',...
    'trialCorrectness' ,'trialResponse','trialRT','propCorrect','stairThresh');

%%
rows = 2;
cols = ceil(length(subFolders)/rows);
for iSub = 1:length(subFolders)
    subplot(rows,cols,iSub)
    for rwd=1:2
        plot(meanPupil{iSub,rwd})
        hold all
    end
    title(expName{iSub,1})
end
figure
for iSub = 1:length(subFolders)
    minlength = min(length(meanPupil{iSub,1}), length(meanPupil{iSub,2}));
    diffPupil{iSub} = meanPupil{iSub,1}(1:minlength) - meanPupil{iSub,2}(1:minlength);
    subplot(rows,cols,iSub)
    plot(diffPupil{iSub});
    hold on
    plot(zeros(1,length(diffPupil{iSub})),'k')
    title(expName{iSub,1})
end
toc