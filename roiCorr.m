close all
clear all
tic
onlyCorrect=0;%1=correct,2=incorrect,0=all
toZscore=0;%0 or 1
regressGlobalMean = 0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
rwdString = {'high','low'};
plotColors = {[1 0 0], [0 0 1], [0 1 0]};
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
globalMeanString = '';
if regressGlobalMean
    globalMeanString = '_globalRegressed';
end

ConcatProjStr = '';
if ConcatProj
    ConcatProjStr = 'ConcatProj';
end
%%
load([dataFolder 'rwdTC_concat' onlyCorrectString zScoreString globalMeanString ConcatProjStr '.mat'], 'concatInfo',  'subResponse', 'roiMeanTseries', ...
    'roiTC', 'allTrials', ...
    'subFolders', 'roiNames','subTrialResponse','trialCorrectness', 'trialResponse', 'trialRT', 'propCorrect',...
    'expName','stairThresh','eccen','ang','areas','trialLength',...
    'subMeanCorrectness', 'subMeanRT','subMedianRT','subMeanThresh',...
    'subMeanRunTC','subStdRunTC','subStd','subRoiRuns',...
    'globalMean','regressBetasGlobal','runRwd',...
    'subRoiRuns','runMeanFFT',...
    'allVoxTrialResponse','allVoxTaskPhase','allVoxTaskAmp','allVoxTaskCo',...
    'voxTrials','voxGoodTrials','meanVoxTrial');
ifig=0;
numSubs=length(subFolders);
goodSubs = [1:3 5:length(subFolders)]; %excluding subject 22
for iSub = 1:length(goodSubs)
    for iRoi=1:length(roiNames)
        for rwd=1:2

            meanVoxTrial{goodSubs(iSub),iRoi,rwd};%vox,t
            stdAmp{iSub,iRoi}(rwd,:) = squeeze(std(meanVoxTrial{goodSubs(iSub),iRoi,rwd},0,2));
            temp = fft(meanVoxTrial{goodSubs(iSub),iRoi,rwd},[],2);
            ph{iSub,iRoi}(rwd,:) = angle(temp(:,2));
            amp{iSub,iRoi}(rwd,:) = abs(temp(:,2));
            
            voxGoodTrials{goodSubs(iSub),iRoi,rwd}; %vox,t,trial
            stdAmpTrial{iSub,iRoi,rwd} = squeeze(std(voxGoodTrials{goodSubs(iSub),iRoi,rwd},0,2));%vox,trial
            temp = fft(voxGoodTrials{goodSubs(iSub),iRoi,rwd},[],2);
            phTrial{iSub,iRoi,rwd} = squeeze(angle(temp(:,2,:)));
            ampTrial{iSub,iRoi,rwd} = squeeze(abs(temp(:,2,:)));
            
            phVar{iSub,iRoi}(rwd,:) = squeeze(circ_std(phTrial{iSub,iRoi,rwd},[],[],2));
            stdAmpVar{iSub,iRoi}(rwd,:) = squeeze(nanstd(stdAmpTrial{iSub,iRoi,rwd},0,2));
            ampVar{iSub,iRoi}(rwd,:) = squeeze(nanstd(ampTrial{iSub,iRoi,rwd},0,2));
            
            
            %ROI average
            reshapedTrials = subTrialResponse{goodSubs(iSub),iRoi,rwd};
            %trial-by-trial FFT
            subTimepointStd(iSub,iRoi,rwd,:) = std(reshapedTrials,0,2);
            temp = fft(reshapedTrials);
            roiFftAmp = abs(temp(2,:));
            roiFftPh = angle(temp(2,:));
            trialRoiFftPhVec{iSub,rwd}(iRoi,:) = roiFftPh;
            trialRoiFftAmpVec{iSub,rwd}(iRoi,:) = roiFftAmp;
            trialRoiStdVec{iSub,rwd}(iRoi,:) = std(reshapedTrials,0,1);
            
        end
                    
        stdAmpDiff{iSub,iRoi} = stdAmp{iSub,iRoi}(1,:) - stdAmp{iSub,iRoi}(2,:);
        phDiff{iSub,iRoi} = circ_dist(ph{iSub,iRoi}(1,:),ph{iSub,iRoi}(2,:));
        ampDiff{iSub,iRoi} = amp{iSub,iRoi}(1,:) - amp{iSub,iRoi}(2,:);
        
        stdAmpVarDiff{iSub,iRoi} = stdAmpVar{iSub,iRoi}(1,:) - stdAmpVar{iSub,iRoi}(2,:);
        phVarDiff{iSub,iRoi} = phVar{iSub,iRoi}(1,:) - phVar{iSub,iRoi}(2,:);        
        ampVarDiff{iSub,iRoi} = ampVar{iSub,iRoi}(1,:) - ampVar{iSub,iRoi}(2,:);      
       
        concatForCorr = [stdAmpDiff{iSub,iRoi};  phDiff{iSub,iRoi}; ampDiff{iSub,iRoi}; ...
            stdAmpVarDiff{iSub,iRoi}; phVarDiff{iSub,iRoi}; ampVarDiff{iSub,iRoi}];
        concatForCorr = concatForCorr(:,~isnan(phDiff{iSub,iRoi}));
        
        [pearsonsCorr(iSub,iRoi,:,:), pearsonsPval(iSub,iRoi,:,:)] = corr(concatForCorr');
    
    end    
end




%%
labels = {'\Deltastd','\Deltaph','\Deltaamp','\Deltavar(std)','\Deltavar(ph)','\Deltavar(amp)'};
ifig=ifig+1; figure(ifig); clf
cols=2;
rows=ceil(length(roiNames)/cols);

meanPearsonsCorr= squeeze(nanmean(pearsonsCorr));%mean over subjects
for iRoi=1:length(roiNames)
   subplot(rows,cols,iRoi)
   imagesc(squeeze(meanPearsonsCorr(iRoi,:,:)));
   xticklabels(labels);
   yticklabels(labels);
   title(roiNames{iRoi});
   axis image
   caxis([-1 1]);
end

% figure(2); clf
% meanPearsonsPval= squeeze(nanmean(pearsonsPval));%mean over subjects
% for iRoi=1:length(roiNames)
%    subplot(rows,cols,iRoi)
%    imagesc(squeeze(meanPearsonsPval(iRoi,:,:)));
%    xticklabels(labels);
%    yticklabels(labels);
%    title(roiNames{iRoi});
% end

%% plot polar histogram of task related phase
ifig=ifig+1; figure(ifig); clf
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

%% Correlate singe trial phases between ROIs , averaged across subjects
%correlate phases across trials between ROIs, within subject & rwd
for iSub = 1:length(goodSubs)
    for rwd=1:2
        %         [phRoiCorr(iSub,rwd,:,:) phRoiPval(iSub,rwd,:,:)] = corr(trialRoiFftPhVec{iSub,rwd}');
        for iRoi1=1:length(roiNames)
            for iRoi2=1:length(roiNames)
                phRoiCorr(iSub,rwd,iRoi1,iRoi2) = circ_corrcc(trialRoiFftPhVec{iSub,rwd}(iRoi1,:),trialRoiFftPhVec{iSub,rwd}(iRoi2,:));
            end
        end
        ampRoiCorr(iSub,rwd,:,:) = corr(trialRoiFftAmpVec{iSub,rwd}');
        stdRoiCorr(iSub,rwd,:,:) = corr(trialRoiStdVec{iSub,rwd}');
    end
end

ifig=ifig+1; figure(ifig); clf
rows=3;
cols=2;
for rwd=1:2
    subplot(rows,cols,rwd);
    imagesc(squeeze(mean(phRoiCorr(:,rwd,:,:))));
    title(['single trial phase | ' rwdString{rwd}]);
    
    subplot(rows,cols,cols+rwd);
    imagesc(squeeze(mean(ampRoiCorr(:,rwd,:,:))));
    title(['single trial FFT amp | ' rwdString{rwd}]);
    
    subplot(rows,cols,2*cols+rwd);
    imagesc(squeeze(mean(stdRoiCorr(:,rwd,:,:))));
    title(['single trial STD | ' rwdString{rwd}]);
    
    for row=1:rows
         subplot(rows,cols,(row-1)*cols+rwd);
            axis square
        caxis([-1 1]);
    end
end

%% Correlate phase of mean response across subjects, between ROIs, within rwd
for iSub = 1:length(goodSubs)
    for iRoi=1:length(roiNames)
        for rwd = 1:2
            subMeanTC(iSub,iRoi,rwd,:) = mean(subTrialResponse{goodSubs(iSub),iRoi,rwd},2);
            temp = fft(squeeze(subMeanTC(iSub,iRoi,rwd,:)));
            subMeanFftAmp(iSub,iRoi,rwd) = abs(temp(2));
            subMeanFftPh(iSub,iRoi,rwd) = angle(temp(2));
        end
    end
end
ifig=ifig+1; figure(ifig); clf
rows=1;
cols=2;
clear roiPhaseCorr
markersize  = 80;
scatterCmap=cool;
% colormap(scatterCmap);
l = size(scatterCmap,1);
for rwd=1:2
    for iRoi1=1:length(roiNames)
        for iRoi2=1:length(roiNames)
            [roiPhaseCorr(rwd,iRoi1,iRoi2), roiPhasePval(rwd,iRoi1,iRoi2)]= circ_corrcc(subMeanFftPh(:,iRoi1,rwd),subMeanFftPh(:,iRoi2,rwd));
        end
        roiPhaseStd(iRoi1,rwd) = circ_std(subMeanFftPh(:,iRoi1,rwd));
        roiColor(iRoi1,:) = scatterCmap(1+floor(iRoi1*(l-1)/(length(roiNames))),:);
    end
    subplot(rows,cols,rwd);
    imagesc(squeeze(roiPhaseCorr(rwd,:,:)));
    axis square
    title(['mean phase per subject | ' rwdString{rwd}]);
    caxis([-1 1]);
end



ifig=ifig+1; figure(ifig); clf
scatter(roiPhaseStd(:,1), roiPhaseStd(:,2),markersize,roiColor,'filled');
% plot(roiPhaseStd(:,1), roiPhaseStd(:,2),'.','markersize',20);
hold on
plot(0:0.1:1.2,0:0.1:1.2,'k-');
colormap(scatterCmap);
colorbar
axis square
xlabel('high')
ylabel('low')
title('phase variability across subjects');

%% Correlate mean task-related response between ROIs, averaged across subjects
for iSub = 1:length(goodSubs)
%     for iRoi=1:length(roiNames)
        for rwd = 1:2
            responseCorr(iSub,rwd,:,:) = corr(squeeze(subMeanTC(iSub,:,rwd,:))');
        end
%     end
end
meanResponseCorr= squeeze(nanmean(responseCorr));%mean over subjects
ifig=ifig+1; figure(ifig); clf
for rwd=1:2
    subplot(rows,cols,rwd)
    imagesc(squeeze(meanResponseCorr(rwd,:,:)));
    title(['mean response timecourse | ' rwdString{rwd}]);
    caxis([-1 1]);
    axis square
end


%%
for iSub = 1:length(goodSubs)
    for iRoi=1:length(roiNames)/2
        subBilateral(iSub,iRoi,:) =squeeze(mean(mean(subMeanTC(iSub,(iRoi-1)*2+1:iRoi*2,:,:),3),2));%mean over both hemispheres, and both rwd
        f(iSub,iRoi,:)=fft(subBilateral(iSub,iRoi,:));
    end
end
subReal = real(f(:,:,2));
subImag = imag(f(:,:,2));

subPh = angle(f(:,:,2));
subAmp = abs(f(:,:,2));

%%
i=i+1; figure(i); clf
titleStr = {'mean','H','L'};
% plot(subReal, subImag, 'k.');
rows=1;cols=3;
nSTD = tinv(0.99,length(subFolders)-1)/sqrt(length(subFolders));
nSTD = 2.5;

twoRois = [1 2; 1 3; 2 3];
rows=size(twoRois,1);
cols=3;

for i=1:size(twoRois,1)
    for iRoi=1:2
        subplot(rows,cols,(i-1)*cols+iRoi)
        scatter(subReal(:,twoRois(i,iRoi)),  subImag(:,twoRois(i,iRoi)), 30, plotColors{twoRois(i,iRoi)},'filled');
        hline(0)
        vline(0)
    end
    subplot(rows,cols,(i-1)*cols+3)
    for iRoi=1:2
        scatter(subReal(:,twoRois(i,iRoi)),  subImag(:,twoRois(i,iRoi)), 30, plotColors{twoRois(i,iRoi)},'filled');
        hold on
    end
    for iSub=1:numSubs
        line( [subReal(:,twoRois(i,1)) subReal(:,twoRois(i,2))]', [subImag(:,twoRois(i,1)) subImag(:,twoRois(i,2))]','color','k');
    end
    
end

for i=1:rows*cols
    subplot(rows,cols,i)
    hline(0)
    vline(0)
end
%%
for iRoi1=1:length(roiNames)/2
    for iRoi2=1:length(roiNames)/2
        [nullSubCorr(iRoi1,iRoi2), nullSubPval(iRoi1,iRoi2)] = circ_corrcc(subPh(:,iRoi1),subPh(:,iRoi2));
    end
end