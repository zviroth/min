onlyCorrect=0;
trialLength=10;
clear dataFolder
dataFolder{1} = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects/';
origFolder = pwd;
subdirs=[];
for i=1:length(dataFolder)
    cd(dataFolder{i});
    temp = dir('s00*');
    if length(temp)==0
        temp = dir('00*');
    end
    subdirs = [subdirs; temp];
end

mrQuit;
if exist('v')
    deleteView(v);
end
cmap = twoCondCmap(64);
cmap = flipud(cmap);
for iSub = 1:length(subdirs)
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
    for iScan = 1:nScans%2 concatenations, 1 for each reward type
        s = viewGet(v, 'stimfile', iScan);
        rwdType = s{1}.stimulus.rewardVal;
        if strcmp(rwdType, 'H')
            rwdTypeNum = 1;
        elseif strcmp(rwdType, 'L')
            rwdTypeNum = 2;
        end
        
        tc{rwdTypeNum} = loadTSeries(v,iScan);
        
        concatInfo{iSub,rwdTypeNum} = viewGet(v, 'concatInfo', iScan);
        
        for r=1:length(s)
            %if subject stopped responding at a certain point we
            %need to fill in the last incorrect trials:
            trialCorrectness{iSub,rwdTypeNum}(r,:) = [s{r}.stimulus.trialCorrectness zeros(1,17-size(s{r}.stimulus.trialCorrectness,2))];
            trialResponse{iSub,rwdTypeNum,r} = s{r}.stimulus.trialResponse;
            trialRT{iSub,rwdTypeNum,r} = s{r}.stimulus.trialRT;
            propCorrect{iSub,rwdTypeNum}(r) = s{r}.stimulus.percentCorrect;
            
        end
        expName{iSub,rwdTypeNum} = s{iScan}.task{1}.taskFilename;
    end
    for rwd=1:2
        [x, y, z, t] = size(tc{rwd});
        voxNum = x*y*z;
        %reshape
        tcReshape{rwd} = reshape(tc{rwd}(:), x,y,z, trialLength, t/trialLength);
        trialResponse{rwd} = mean(tcReshape{rwd},5);
        responseAmp{rwd} = std(trialResponse{rwd},0,4,'omitnan');
        %             tcReshape{rwd} = reshape(tc{rwd}(:), voxNum, trialLength, t/trialLength);
        %             trialResponse{rwd} = mean(tcReshape{rwd},3);
        
        
%         average per run - ALL TRIALS AVERAGED, NOT ONLY CORRECT TRIALS
        tcReshapeRun{rwd} = reshape(tcReshape{rwd},x,y,z,10,15,[]);
        runResponse{rwd} = squeeze(mean(tcReshapeRun{rwd},5));
            
        runFFT{rwd} = fft(permute(runResponse{rwd},[4 1 2 3 5]));
        runFFTamp{rwd} = squeeze(abs(runFFT{rwd}(2,:,:,:,:)));
        runFFTphase{rwd} = squeeze(angle(runFFT{rwd}(2,:,:,:,:)));
        [s s0] = circ_std(runFFTphase{rwd},[],[],4);
        stdRunFFTphase{rwd} = s0;
        meanRunFFTamp{rwd} = mean(runFFTamp{rwd},4);
        stdRunFFTamp{rwd} = std(runFFTamp{rwd},0,4);
    end
    
    
%     timeStdDiff = stdRunFFTphase{1}-stdRunFFTphase{2};
    timeStdDiffInd = (stdRunFFTphase{1}-stdRunFFTphase{2})./(stdRunFFTphase{1}+stdRunFFTphase{2});
%     fftAmpDiff = meanRunFFTamp{1} - meanRunFFTamp{2};
%     fftAmpStdDiff = stdRunFFTamp{1} - stdRunFFTamp{2};
    fftAmpStdDiffInd =  (stdRunFFTamp{1} - stdRunFFTamp{2}) ./ (stdRunFFTamp{1} + stdRunFFTamp{2});
%     responseDiff = responseAmp{1}-responseAmp{2};
    respDiffInd = (responseAmp{1}-responseAmp{2})./(responseAmp{1}+responseAmp{2});
    mrSetPref('overwritePolicy','Overwrite');
    overlayScan = 1;
    groupNum = viewGet(v, 'curGroup');
%     [v rwdOverlay] = mrDispOverlay(responseDiff,overlayScan,groupNum,v,'overlayName=amplitudeDiff','saveName=rwdOverlay', 'cmap', parula);
    [v rwdOverlay] = mrDispOverlay(respDiffInd,overlayScan,groupNum,v,'overlayName=amplitudeDiff','saveName=rwdOverlay', 'cmap', cmap,'colormapType=setRangeToMaxAroundZero');
    [v rwdOverlay] = mrDispOverlay(timeStdDiffInd,overlayScan,groupNum,v,'overlayName=timingStdDiff','saveName=rwdOverlay', 'cmap', cmap,'colormapType=setRangeToMaxAroundZero');
    [v rwdOverlay] = mrDispOverlay(fftAmpStdDiffInd,overlayScan,groupNum,v,'overlayName=fftAmpStdDiff','saveName=rwdOverlay', 'cmap', cmap,'colormapType=setRangeToMaxAroundZero');
    
    deleteView(v);
    cd('..');
end

cd(origFolder);

