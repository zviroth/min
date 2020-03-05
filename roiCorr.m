close all
clear all
tic
onlyCorrect=0;%1=correct,2=incorrect,0=all trials with response, 4=all trials.
toZscore=1;%0 or 1
regressGlobalMean = 0;
ConcatProj = 1;
curFolder = pwd;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
rwdString = {'high','low'};
plotColors = {[1 0 0], [0 0 1], [0 1 0], [1 1 0]};
%%
onlyCorrectString = '';
if onlyCorrect==1
    onlyCorrectString = '_correct';
elseif onlyCorrect==2
    onlyCorrectString = '_incorrect';
elseif onlyCorrect==0
    onlyCorrectString = '_validresponse';
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
    'voxTrials','voxGoodTrials','meanVoxTrial',...
    'maxRT');
roiLabels = {'benson','localizer','DMN','cerebellum'};
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
            
            temp= squeeze(std(voxGoodTrials{goodSubs(iSub),iRoi,rwd},0,3));%vox,timepoint
            subVoxVar{iSub,iRoi}(rwd,:) = mean(temp,2);
            
            %ROI average
            reshapedTrials = subTrialResponse{goodSubs(iSub),iRoi,rwd};
            %trial-by-trial FFT
            subTimepointStd(iSub,iRoi,rwd,:) = std(reshapedTrials,0,2);%timepoint variability timecourse (=std across trials)
            temp = fft(reshapedTrials);
            roiFftAmp = abs(temp(2,:));
            roiFftPh = angle(temp(2,:));
            trialRoiFftPhVec{iSub,rwd}(iRoi,:) = roiFftPh;
            trialRoiFftAmpVec{iSub,rwd}(iRoi,:) = roiFftAmp;
            trialRoiStdVec{iSub,rwd}(iRoi,:) = std(reshapedTrials,0,1);
            
            subMeanStd(iSub,iRoi,rwd) = std(mean(reshapedTrials,2));%std of mean response
            subMeanVar(iSub,iRoi,rwd) = mean(squeeze(subTimepointStd(iSub,iRoi,rwd,:)));%mean timepoint variability
            
        end
                    
        stdAmpDiff{iSub,iRoi} = stdAmp{iSub,iRoi}(1,:) - stdAmp{iSub,iRoi}(2,:);
        phDiff{iSub,iRoi} = circ_dist(ph{iSub,iRoi}(1,:),ph{iSub,iRoi}(2,:));
        ampDiff{iSub,iRoi} = amp{iSub,iRoi}(1,:) - amp{iSub,iRoi}(2,:);
        
        stdAmpVarDiff{iSub,iRoi} = stdAmpVar{iSub,iRoi}(1,:) - stdAmpVar{iSub,iRoi}(2,:);
        phVarDiff{iSub,iRoi} = phVar{iSub,iRoi}(1,:) - phVar{iSub,iRoi}(2,:);        
        ampVarDiff{iSub,iRoi} = ampVar{iSub,iRoi}(1,:) - ampVar{iSub,iRoi}(2,:);      
        timepointVarDiff{iSub,iRoi} = subVoxVar{iSub,iRoi}(1,:) - subVoxVar{iSub,iRoi}(2,:);      
%         concatForCorr = [stdAmpDiff{iSub,iRoi};  phDiff{iSub,iRoi}; ampDiff{iSub,iRoi}; ...
%             stdAmpVarDiff{iSub,iRoi}; phVarDiff{iSub,iRoi}; ampVarDiff{iSub,iRoi}];
%         concatForCorr = [stdAmpDiff{iSub,iRoi};   ampDiff{iSub,iRoi}; ...
%             phDiff{iSub,iRoi}; stdAmpVarDiff{iSub,iRoi}; ampVarDiff{iSub,iRoi}; phVarDiff{iSub,iRoi}];
        concatForCorr = [mean(stdAmp{iSub,iRoi}); mean(amp{iSub,iRoi}); circ_mean(ph{iSub,iRoi}); ...
            stdAmpDiff{iSub,iRoi};   ampDiff{iSub,iRoi}; phDiff{iSub,iRoi}; ...
            mean(stdAmpVar{iSub,iRoi}); mean(phVar{iSub,iRoi}); mean(ampVar{iSub,iRoi}); mean(subVoxVar{iSub,iRoi}); ...
             stdAmpVarDiff{iSub,iRoi}; ampVarDiff{iSub,iRoi}; phVarDiff{iSub,iRoi}; timepointVarDiff{iSub,iRoi}];
        
        
        concatForCorr = concatForCorr(:,~isnan(phDiff{iSub,iRoi}));
        
        [pearsonsCorr(iSub,iRoi,:,:), pearsonsPval(iSub,iRoi,:,:)] = corr(concatForCorr');
    
    end    
end



%%
% labels = {'\Deltastd','\Deltaph','\Deltaamp','\Deltavar(std)','\Deltavar(ph)','\Deltavar(amp)'};
labels = {'std', 'amp','ph', '\Deltastd','\Deltaamp','\Deltaph','std var', 'amp var', 'ph var', 'time var','\Deltavar(std)','\Deltavar(amp)','\Deltavar(ph)','\Deltavar(time)'};
ifig=ifig+1; figure(ifig); clf
cols=2;
rows=ceil(length(roiNames)/cols);

meanPearsonsCorr= squeeze(nanmean(pearsonsCorr));%mean over subjects
for iRoi=1:length(roiNames)
   subplot(rows,cols,iRoi)
   imagesc(squeeze(meanPearsonsCorr(iRoi,:,:)));
   xticklabels(labels(1:3:end));
   xticks(1:3:length(labels));
   yticklabels(labels);
   yticks(1:length(labels));
   if iRoi==1
       title({'mean corr across vox'; roiNames{iRoi}});
   else
   title(roiNames{iRoi});
   end
   axis image
   caxis([-1 1]);
%    title(['mean corr between voxels']);
end
 set(gcf,'position',[50 50 500 	800]);
% figure(2); clf
% meanPearsonsPval= squeeze(nanmean(pearsonsPval));%mean over subjects
% for iRoi=1:length(roiNames)
%    subplot(rows,cols,iRoi)
%    imagesc(squeeze(meanPearsonsPval(iRoi,:,:)));
%    xticklabels(labels);
%    yticklabels(labels);
%    title(roiNames{iRoi});
% end

% %% plot polar histogram of task related phase
% ifig=ifig+1; figure(ifig); clf
% rows=2;
% cols=numSubs;
% nbins=20;
% nVox=1000;
% thresh=0.015;
% for iSub=1:numSubs
% %     bestVox = sort(allVoxTaskCo{iSub,1}+allVoxTaskCo{iSub,2},'descend');
%     for rwd=1:2
%         subplot(rows,cols,iSub+(rwd-1)*cols)
% %         polarhistogram(allVoxTaskPhase{iSub,rwd}(bestVox(1:nVox)),nbins);
%         polarhistogram(allVoxTaskPhase{iSub,rwd}(allVoxTaskCo{iSub,rwd}>thresh),nbins);
%         set(gca,'RTickLabel',[]);
%         set(gca,'ThetaTickLabel',[]);
%     end
% end

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
for isubplot=1:cols*rows
    subplot(rows,cols,isubplot)
    yticklabels(roiNames);
    yticks(1:length(roiNames));
    xticklabels(roiLabels);
%     xticks(1:length(roiNames));
end
set(gcf,'position',[120 120 800 800]);
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

    yticklabels(roiNames);
    yticks(1:length(roiNames));
    xticklabels(roiLabels);
%     xticks(1:length(roiNames));

end
set(gcf,'position',[220 220 800 400]);


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
rows=1;
cols=2;
for rwd=1:2
    subplot(rows,cols,rwd)
    imagesc(squeeze(meanResponseCorr(rwd,:,:)));
    title(['mean response timecourse | ' rwdString{rwd}]);
    caxis([-1 1]);
    axis square
    yticklabels(roiNames)
    xticklabels(roiLabels);
end
set(gcf,'position',[120 120 800 400]);

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
ifig=ifig+1; figure(ifig); clf
titleStr = {'mean','H','L'};
% plot(subReal, subImag, 'k.');
rows=1;cols=3;
nSTD = tinv(0.99,length(subFolders)-1)/sqrt(length(subFolders));
nSTD = 2.5;

twoRois = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

rows=size(twoRois,1);
cols=3;

for i=1:size(twoRois,1)
    for iRoi=1:2
        subplot(rows,cols,(i-1)*cols+iRoi)
        scatter(subReal(:,twoRois(i,iRoi)),  subImag(:,twoRois(i,iRoi)), 30, plotColors{twoRois(i,iRoi)},'filled');
        hline(0)
        vline(0)
        title(roiLabels{twoRois(i,iRoi)});
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
    axis square
end
 set(gcf,'position',[50 50 500 	800]);
%%
for iRoi1=1:length(roiNames)/2
    for iRoi2=1:length(roiNames)/2
        [nullSubCorr(iRoi1,iRoi2), nullSubPval(iRoi1,iRoi2)] = circ_corrcc(subPh(:,iRoi1),subPh(:,iRoi2));
    end
end

%% ECCENTRICITY BINS - ONLY FOR BENSON ROIs
%bins by log eccentricity
eccMin = 0.2;
eccMax = 70;
nbins = 12;
% binBorders = (linspace(eccMin^3,eccMax^3,nbins+1))^(1/3)
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
% binBorders = [0.2    0.9    1.55    2.45    3.7    5    6.95   10.1   18   35   70.0000];
% binBorders = [0.2 3 70];


nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end

clear subBinTrialResponse subBinResponse binVoxels numBinVoxels meanSubVoxCorr meanVoxCorr
for iSub = 1:length(goodSubs)%length(subdirs)
    for iRoi=1:2%length(roiNames)
        for rwd=1:2
            temp = trialCorrectness{goodSubs(iSub),rwd}(:,2:end-1);
            trialCorrectnessVec = temp(:);
            temp = trialResponse{goodSubs(iSub),rwd}(:,2:end-1);
            trialResponseVec = temp(:);
            temp = trialRT{goodSubs(iSub),rwd}(:,2:end-1);
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
                
            for ibin=1:nbins
                binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
%                 binVoxels = binVoxels & areas{iSub,iRoi}==1;%ONLY V1
                binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries = nanmean(roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:));%mean timecourse across voxels
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = reshape(binMeanTseries, trialLength, length(binMeanTseries)/trialLength);
                
                
                subBinTrialResponse{iSub,iRoi,ibin,rwd} = subBinTrialResponse{iSub,iRoi,ibin,rwd}(:,goodTrials);
                
                
                reshapedTrials = reshape(subBinTrialResponse{iSub,iRoi,ibin,rwd},trialLength,[]);
                if iSub==1
                    allBinTrials{iRoi,ibin,rwd} = reshapedTrials;
                else
                    allBinTrials{iRoi,ibin,rwd} = [allBinTrials{iRoi,ibin,rwd} reshapedTrials];
                end
                %average per subject
                subBinResponse(iSub,iRoi,ibin,rwd,:) = mean(subBinTrialResponse{iSub,iRoi,ibin,rwd},2);%mean response across trials
                %trial-by-trial variability per subject
                subBinVar(iSub,iRoi,ibin,rwd) = mean(std(subBinTrialResponse{iSub,iRoi,ibin,rwd},0,2));%mean timepoint variability
                
                temp = fft(subBinTrialResponse{iSub,iRoi,ibin,rwd});
                roiBinFftAmp{iSub,iRoi,rwd}(ibin,:) = abs(temp(2,:));
                roiBinFftPh{iSub,iRoi,rwd}(ibin,:) = angle(temp(2,:));
%                 subBinFftAmpVar(iSub,iRoi,ibin,rwd) = std(roiBinFftAmp{iSub,iRoi,rwd}(ibin,:));
%                 subBinFftPhVar(iSub,iRoi,ibin,rwd) = circ_std(roiBinFftPh{iSub,iRoi,rwd}(ibin,:)');
                
                subBinStd{iSub,iRoi,rwd}(ibin,:) = std(subBinTrialResponse{iSub,iRoi,ibin,rwd});%response std per trial
%                 singleTrialStd = std(subBinTrialResponse{iSub,iRoi,ibin,rwd});
%                 subBinStdVar(iSub,iRoi,ibin,rwd) = std(subBinStd{iSub,iRoi,rwd}(ibin,:));
%                 subBinStdMean(iSub,iRoi,ibin,rwd) = mean(subBinStd{iSub,iRoi,rwd}(ibin,:));
                temp = fft(squeeze(subBinResponse(iSub,iRoi,ibin,rwd,:)));
                subBinMeanStd(iSub,iRoi,ibin,rwd) = std(squeeze(subBinResponse(iSub,iRoi,ibin,rwd,:)));
                subBinStdVar(iSub,iRoi,ibin,rwd) = std(subBinStd{iSub,iRoi,rwd}(ibin,:));
                subBinMeanPh(iSub,iRoi,ibin,rwd) = angle(temp(2));
                subBinPhVar(iSub,iRoi,ibin,rwd) = circ_std(roiBinFftPh{iSub,iRoi,rwd}(ibin,:)');
                subBinMeanAmp(iSub,iRoi,ibin,rwd) = abs(temp(2));
                subBinAmpVar(iSub,iRoi,ibin,rwd) = std(roiBinFftAmp{iSub,iRoi,rwd}(ibin,:));

            end
            %correlation between bins
            for ibin1=1:nbins
                for ibin2 = 1:nbins
%                     subBinCorrPh(iSub,iRoi, rwd, ibin1, ibin2) = circ_corrcc(roiBinFftPh{iSub,iRoi,rwd}(ibin1,:), roiBinFftPh{iSub,iRoi,rwd}(ibin2,:));
                    subBinCorrPh(iSub,iRoi, rwd, ibin1, ibin2) = corr(roiBinFftPh{iSub,iRoi,rwd}(ibin1,:)', roiBinFftPh{iSub,iRoi,rwd}(ibin2,:)');
                end
            end
            subBinCorrAmp(iSub,iRoi, rwd,:,:) = corr(roiBinFftAmp{iSub,iRoi,rwd}');
        end
    end

    %add other ROIs
    for rwd=1:2
        for iRoi=3:length(roiNames)
            ihemi = 2-mod(iRoi,2);%1==left, 2==right
            iAddRoi = nbins+ceil((iRoi-2)/2);
            %value per trial for each subject/roi
            roiBinFftPh{iSub,ihemi,rwd}(iAddRoi,:) = trialRoiFftPhVec{iSub,rwd}(iRoi,:);
            roiBinFftAmp{iSub,ihemi,rwd}(iAddRoi,:) = trialRoiFftAmpVec{iSub,rwd}(iRoi,:);
            subBinStd{iSub,ihemi,rwd}(iAddRoi,:) = trialRoiStdVec{iSub,rwd}(iRoi,:);
            %value per subject/roi
            subBinMeanStd(iSub,ihemi,iAddRoi,rwd) = subMeanStd(iSub,iRoi,rwd);
            subBinMeanPh(iSub,ihemi,iAddRoi,rwd) = subMeanFftPh(iSub,iRoi,rwd);
            subBinMeanAmp(iSub,ihemi,iAddRoi,rwd) = subMeanFftAmp(iSub,iRoi,rwd);
            subBinStdVar(iSub,ihemi,iAddRoi,rwd) = subMeanVar(iSub,iRoi,rwd);%should be same as std(trialRoiStdVec{iSub,rwd}(iRoi,:))
            subBinPhVar(iSub,ihemi,iAddRoi,rwd) = circ_std(trialRoiFftPhVec{iSub,rwd}(iRoi,:)');
            subBinAmpVar(iSub,ihemi,iAddRoi,rwd) = std(trialRoiFftAmpVec{iSub,rwd}(iRoi,:));

%             subBinVar(iSub,iRoi,end+1,rwd) = mean(subTimepointStd(iSub,iRoi,rwd,:));%mean timepoint variability
            
        end
        %     trialRoiFftPhVec{iSub,rwd}(iRoi1,:)
    end
    
    %correlate between hemispheres
    for rwd=1:2
        for ibin1=1:size(roiBinFftPh{iSub,1,rwd},1)
            for ibin2 = 1:size(roiBinFftPh{iSub,2,rwd},1)
                subHemiBinCircCorrPh(iSub, rwd, ibin1, ibin2) = circ_corrcc(roiBinFftPh{iSub,1,rwd}(ibin1,:), roiBinFftPh{iSub,2,rwd}(ibin2,:));
                subHemiBinCorrPh(iSub, rwd, ibin1, ibin2) = corr(roiBinFftPh{iSub,1,rwd}(ibin1,:)', roiBinFftPh{iSub,2,rwd}(ibin2,:)');

                hemiBinCorrPh(rwd,ibin1,ibin2) = circ_corrcc(subBinMeanPh(:,1,ibin1,rwd),subBinMeanPh(:,2,ibin2,rwd));
            end
        end
        subHemiBinCorrAmp(iSub, rwd,:,:) = corr(roiBinFftAmp{iSub,1,rwd}', roiBinFftAmp{iSub,2,rwd}');
        subHemiBinCorrStd(iSub, rwd,:,:) = corr(subBinStd{iSub,1,rwd}', subBinStd{iSub,2,rwd}');

    end
end
for rwd=1:2
   hemiBinCorrAmp(rwd,:,:) = corr(squeeze(subBinMeanAmp(:,1,:,rwd)),squeeze(subBinMeanAmp(:,2,:,rwd))); 
   hemiBinCorrStd(rwd,:,:) = corr(squeeze(subBinMeanStd(:,1,:,rwd)),squeeze(subBinMeanStd(:,2,:,rwd))); 
   hemiBinCorrStdVar(rwd,:,:) = corr(squeeze(subBinStdVar(:,1,:,rwd)),squeeze(subBinStdVar(:,2,:,rwd)));
   hemiBinCorrPhVar(rwd,:,:) = corr(squeeze(subBinPhVar(:,1,:,rwd)),squeeze(subBinPhVar(:,2,:,rwd))); 
   hemiBinCorrAmpVar(rwd,:,:) = corr(squeeze(subBinAmpVar(:,1,:,rwd)),squeeze(subBinAmpVar(:,2,:,rwd)));
end
%% correlations of trial phases between bins, separately for each hemisphere
ifig=ifig+1; figure(ifig); clf
rows=2;
cols=2;
for iRoi=1:2%length(roiNames)
    for rwd=1:2
        subplot(rows,cols,iRoi + (rwd-1)*cols);
        imagesc(squeeze(mean(subBinCorrPh(:,iRoi,rwd,:,:))));
        title(['latency: ' roiNames{iRoi} ' ' rwdString{rwd}]);
        axis square
    end
end
%% correlations of trial phases between bins, between hemispheres
ifig=ifig+1; figure(ifig); clf
rows=4;
cols=2;
%linear correlations of trial phases
for rwd=1:2
    subplot(rows,cols,rwd);
    imagesc(squeeze(mean(subHemiBinCorrPh(:,rwd,:,:))));
    title({'corr across trials'; [ 'latency, linear,  ' rwdString{rwd}]});
    xticks;
    xticklabels;
end
%circular correlations of trial phases
for rwd=1:2
    subplot(rows,cols,rwd+cols);
    imagesc(squeeze(mean(subHemiBinCircCorrPh(:,rwd,:,:))));
    title(['latency, circular,  ' rwdString{rwd}]);
end
% correlations of trial FFT amplitudes between bins, between hemispheres
for rwd=1:2
    subplot(rows,cols,rwd+2*cols);
    imagesc(squeeze(mean(subHemiBinCorrAmp(:,rwd,:,:))));
    title(['fft amplitude,  '  rwdString{rwd}]);
end
% correlations of trial STD amplitudes between bins, between hemispheres
for rwd=1:2
    subplot(rows,cols,rwd+3*cols);
    imagesc(squeeze(mean(subHemiBinCorrStd(:,rwd,:,:))));
    title(['std amplitude,  '  rwdString{rwd}]);
end
for isubplot=1:rows*cols
   subplot(rows,cols,isubplot)
      xlabel('right'); ylabel('left');
   caxis([-1 1]);
   axis square
   colormap jet
end

set(gcf,'position',[450 150 400 	900]);
 
 
 %% correlations across subjects, between hemispheres
ifig=ifig+1; figure(ifig); clf
rows=6;
cols=2;
%mean phase
for rwd=1:2
    subplot(rows,cols,rwd);
    imagesc(squeeze(hemiBinCorrPh(rwd,:,:)));
    title({'corr across subjects'; ['latency, circular,  ' rwdString{rwd}]});
end
%mean amp
for rwd=1:2
    subplot(rows,cols,rwd+cols);
    imagesc(squeeze(hemiBinCorrAmp(rwd,:,:)));
    title(['fft amplitude,  ' rwdString{rwd}]);
end
%mean std amplitude
for rwd=1:2
    subplot(rows,cols,rwd+2*cols);
    imagesc(squeeze(hemiBinCorrStd(rwd,:,:)));
    title(['std,  ' rwdString{rwd}]);
end
%phase variability
for rwd=1:2
    subplot(rows,cols,rwd+3*cols);
    imagesc(squeeze(hemiBinCorrPhVar(rwd,:,:)));
    title(['phase variability,  ' rwdString{rwd}]);
end
%fft amplitude variability
for rwd=1:2
    subplot(rows,cols,rwd+4*cols);
    imagesc(squeeze(hemiBinCorrAmpVar(rwd,:,:)));
    title(['amplitude variability,  ' rwdString{rwd}]);
end
%std amplitude variability
for rwd=1:2
    subplot(rows,cols,rwd+5*cols);
    imagesc(squeeze(hemiBinCorrStdVar(rwd,:,:)));
    title(['std variability,  ' rwdString{rwd}]);
end



for isubplot=1:rows*cols
   subplot(rows,cols,isubplot)
   xlabel('right'); ylabel('left');
   xticklabels('');
   xticks('');
   yticklabels('');
   yticks('');
   caxis([-1 1]);
   axis square
   colormap jet
end
set(gcf,'position',[450 150 400 	900]);
%%
ifig=ifig+1; figure(ifig); clf
temp = squeeze(mean(numBinVoxels));
plot(mean(temp));