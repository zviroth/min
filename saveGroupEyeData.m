close all
clear all
tic
ConcatProj = 1;
saveFolder = '/Volumes/MH02086153MACDT-Drobo/allMinBehavioral/';
dataFolder = saveFolder;
samplerate=500;
subFolders = {'s000520180126','s000820171116','s001620171103','s002220171122','s003220180105','s003920180220',...
    's004020180221','s004120180223','s004320180306','s004520180308','s004620180323','s004720180326',...
    's004920180404'};

numSubs = length(subFolders);
curFolder = pwd;
onlyCorrect=0;
cd(dataFolder);
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
clear d concatInfo sumResponse subResponse roiMeanTseries meanResponse roiTC subRunResponse trialCorrectness trialResponse trialRT propCorrect
clear contrast spatFreq nullTrial contrastConcat spatFreqConcat nullTrialConcat designMatContrast designMatFreq nullDeconv deconvMatNull meanDeconvNull meanDeconvNullConstant
clear taskResp taskPrediction meanPupil
clear ecg ecgPeaks ecgPeaksDiff ecgPeaksAmp resp respPeaks respPeaksDiff respPeaksAmp
clear designMatContrastFreq designMatNull taskBaseline taskTimes rwdPupil 
clear sacRun pupilRun startTimes endTimes eyetrackerTime listRwdRuns

%smooth out the microsaccade data
L = 100;
filter = ones(L,1);
% filter = [1:L/2 L/2:-1:1];
filter = filter./sum(filter);
    
for iSub = 1:numSubs
    cd(subFolders{iSub});
    
    stimfiles = dir('*.mat');
    listRwdRuns{iSub,1} = [];
    listRwdRuns{iSub,2} = [];
    for iRun=1:length(stimfiles)
        stimfile=stimfiles(iRun).name;
        allStimfiles{iRun}=load(stimfile);
        
        expName{iSub,iRun} = allStimfiles{iRun}.task{1}.taskFilename;

        rwdType = allStimfiles{iRun}.stimulus.rewardVal;
        if strcmp(rwdType, 'H')
            rwd = 1;
        elseif strcmp(rwdType, 'L')
            rwd = 2;
        else
            disp('wtf');
            keyboard
        end
        concatRwdTypeNum(iSub,iRun) = rwd;
        listRwdRuns{iSub,rwd} = [listRwdRuns{iSub,rwd}; iRun];
        
        
    end
    for rwd=1:2
        clear s
        for r=1:length(listRwdRuns{iSub,rwd})%going through all stimfiles for this reward level
            iRun = listRwdRuns{iSub,rwd}(r);
            s{r} = allStimfiles{iRun};
            trialCorrectness{iSub,rwd}(r,:) = [s{r}.stimulus.trialCorrectness zeros(1,17-size(s{r}.stimulus.trialCorrectness,2))];
            trialResponse{iSub,rwd}(r,:) = [s{r}.stimulus.trialResponse zeros(1,17-size(s{r}.stimulus.trialResponse,2))];%0 means no response
            trialRT{iSub,rwd}(r,:) = [s{r}.stimulus.trialRT NaN(1,17-size(s{r}.stimulus.trialRT,2))];
            propCorrect{iSub,rwd}(r) = s{r}.stimulus.percentCorrect;
            stairThresh{iSub,rwd}(r,:) = [s{r}.stimulus.stair{1}.threshold s{r}.stimulus.stair{2}.threshold];
                
            stimfile = stimfiles(iRun).name;
%             stimfile = s{r}.filename;
            e{r} = myGetTaskEyeTraces(stimfile, 'removeBlink=3');
            startTimes{iSub,rwd,r} = e{r}.complete.startTime;
            endTimes{iSub,rwd,r} = e{r}.complete.endTime;
            eyetrackerTime{iSub,rwd,r} = e{r}.complete.time;
            runEyeSize{rwd}(r,:) = size(e{r}.eye.pupil);
        end
        trialLengthEye(rwd) = max(runEyeSize{rwd}(:,2));
        numTrials(rwd) = sum(runEyeSize{rwd}(:,1));
        sacRwd{iSub,rwd} = zeros(numTrials(rwd), trialLengthEye(rwd));
        trialCounter=0;
        numRuns(iSub,rwd) = length(s);
        for r=1:length(s)
            
            eyePos = [e{r}.complete.xPos' e{r}.complete.yPos'];% time X 2, x and y
            if sum(~isnan(eyePos(:)))~=0
                eyevel{r} = vecvel(eyePos,samplerate, 2);
                sacOutput = microsacc(eyePos, eyevel{r}, 1, 5);
                sacOnset = sacOutput(:,1);               
                sacRun{iSub,rwd,r} = sacOnset;%saccade onset
                sacRunTC = zeros(size(eyePos,1),1);
                sacRunTC(sacOnset) = ones;
                sacRunSmooth{iSub,rwd,r} = conv(sacRunTC, filter, 'same');
                pupilRun{iSub,rwd,r} = e{r}.complete.pupil;%e{r}.eye.pupil is longer because of zero padding trials;

                for itrial=1:length(e{r}.complete.startTime)
                    trialCounter=trialCounter+1;
                    sac{r,itrial} = sacOnset(sacOnset >= e{r}.complete.startTime(itrial) & sacOnset<=e{r}.complete.endTime(itrial));%saccade onsets for this trial
                    sac{r,itrial} = sac{r,itrial} - e{r}.complete.startTime(itrial) + 1;%time relative to trial
                    sacTrial = zeros(1,trialLengthEye(rwd));
                    sacTrial(sac{r,itrial}) = ones;
                    sacRwd{iSub,rwd}(trialCounter,:) = sacTrial;
                    
                    sacRunSmoothTrial{r,itrial} = sacRunSmooth{iSub,rwd,r}(e{r}.complete.startTime(itrial):e{r}.complete.endTime(itrial));
%                     sacTrial = zeros(1,trialLengthEye(rwd));
                    sacTrial = NaN(1,trialLengthEye(rwd));
                    sacTrial(1:length(sacRunSmoothTrial{r,itrial})) = sacRunSmoothTrial{r,itrial};
                    sacRunSmoothTrialZeroFilled{iSub,rwd}(trialCounter,:) = sacTrial;
                end
            end
        end

        binnedSac{iSub,rwd} = sum(sacRwd{iSub,rwd});
        rwdPupil{iSub,rwd} = NaN(numTrials(rwd), trialLengthEye(rwd));
        nullTrials{iSub,rwd} = NaN(numTrials(rwd),1);
        trialCounter=0;
        for r=1:length(s)         
            rwdPupil{iSub,rwd}(trialCounter+1:trialCounter+runEyeSize{rwd}(r,1), 1:runEyeSize{rwd}(r,2)) = e{r}.eye.pupil;
            trialCounter = trialCounter + runEyeSize{rwd}(r,1);
        end
        %plot mean pupil size
        meanPupil{iSub,rwd} = nanmean(rwdPupil{iSub,rwd})';  
        
        smoothSac{iSub,rwd} = nanmean(sacRunSmoothTrialZeroFilled{iSub,rwd});
    end
    
%     keyboard

%     for rwd=1:2
%         smoothSac{iSub,rwd} = conv(binnedSac{iSub,rwd}, filter, 'same');
%     end

    
    %% return to home directory
    cd(dataFolder);
end

%%
save([dataFolder 'eyeData.mat'], 'dataFolder', 'subFolders', 'samplerate',  ...
    'numSubs', 'onlyCorrect', ...
    'sacRwd','binnedSac', 'smoothSac', ...
    'rwdPupil','meanPupil',...
    'sacRun','pupilRun', 'numRuns','startTimes','endTimes','eyetrackerTime',...
    'listRwdRuns', 'concatRwdTypeNum','expName',...
    'sacRunSmooth', 'sacRunSmoothTrialZeroFilled','sacRunSmoothTrial','L');

toc
%%
figure(1); clf;
rows=2;
cols=numSubs;
for iSub=1:numSubs
    subplot(rows,cols,iSub)
    for rwd=1:2
        plot(smoothSac{iSub,rwd})
        hold on
    end
end

for iSub=1:numSubs
    subplot(rows,cols,iSub+cols)
    for rwd=1:2
        plot(binnedSac{iSub,rwd})
        hold on
    end
end