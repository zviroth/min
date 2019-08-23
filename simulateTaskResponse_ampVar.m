close all
clear all
upsampleFactor = 10;
TR = 1.5;
%set canonical HRF
modelParams = struct;
sampleDuration = TR/upsampleFactor;
sampleDelay=sampleDuration/2;
defaultParams=1;
% params that make the variability timecourse similar to temporal jitter
modelParams.x = 12;
modelParams.y = 16;
modelParams.z = 6;

[modelParams hrfModel] = hrfDoubleGamma(modelParams,sampleDuration,sampleDelay,defaultParams);
i=0;
plot(hrfModel);
% keyboard

%%
plotColorsSnr = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2]};
numTrials = 16;
trialLength = 10;
runsPerRwd = 100;
T = numTrials*trialLength;
upT = T*upsampleFactor;
oneOverF(1,1,1,:) = (log10(1:T))*1./(1:T);
oneOverF(1,1,1,2:end) = oneOverF(1,1,1,1:end-1);%we only care about the second component onwards
oneOverF(1,1,1,T/2+1:end) = oneOverF(1,1,1,1+T/2:-1:2);%make it symmetric


i=i+1; figure(i) ;clf
subplot(1,2,1)
plot(hrfModel(1:upsampleFactor:end))
hold on
plot(abs(diff(hrfModel(1:upsampleFactor:end))))
vline(10);

subplot(1,2,2)
plot(squeeze(oneOverF))


%%

% snr = [1 0.95];
% temporalJitter = [0 1.5];%std of jitter, in Fourier phase
npoints = 10;
% nsr = linspace(0.01, 5, npoints);
% snr = 1./nsr;
% snr = linspace(5, 0.01, npoints);
snr = logspace(2, -1, npoints);
% logspace(-2,1,5)
% snr = [3:-0.1:0.1];
% jitter = zeros(size(snr));

jitter = linspace(0, pi/2, npoints);

% jitter = [0 : 0.2: 2*pi];%std of jitter, in Fourier phase
% snr = 1*ones(size(temporalJitter));

% ampJitter = linspace(0, 2, npoints);
%from 0.1 to 0.5 works

ampJitter = logspace(-2,0,npoints);
% ampJitter = 1./linspace(4,0.8,npoints);
% ampJitter = linspace(0.05,0.45,npoints);

for isnr=1:length(snr)
    plotColorsSnr{isnr} = [1 - (isnr-1)/(length(snr)-1), 0, (isnr-1)/(length(snr)-1)];
end
for ijitter=1:length(jitter)
    plotColorsJitter{ijitter} = [1 - (ijitter-1)/(length(jitter)-1), 0, (ijitter-1)/(length(jitter)-1)];
end
for iampJitter=1:length(ampJitter)
    plotColorsAmpJitter{iampJitter} = [1 - (iampJitter-1)/(length(ampJitter)-1), 0, (iampJitter-1)/(length(ampJitter)-1)];
end
% snr = [1 1];
% temporalJitter = [0 1.5];%std of jitter, in Fourier phase



taskAmp = 1;%ones(numRwd,1);


for isnr=1:length(snr)
    for ijitter = 1:length(jitter)
        for iampJitter=1:length(ampJitter)
            for r=1:runsPerRwd
                runTC = zeros(1,upT);
                taskTiming = 1:trialLength*upsampleFactor:upT;%temporal jitter will be added here, after upsampling
                noisyTiming = taskTiming + (jitter(ijitter)*upsampleFactor*trialLength/(2*pi))*randn(size(taskTiming));
                noisyTiming(noisyTiming<1) = 1;
                
                noisyTiming(noisyTiming>T*upsampleFactor) = T;
                
                runTC(ceil(noisyTiming)) = ones;
                runTC(ceil(noisyTiming)) = taskAmp*(1+atan(ampJitter(iampJitter)*randn(numTrials,1)));
                temp = conv(runTC,hrfModel);
                runTC = temp(1:upT);%crop end
                
%                 rwdSignal(isnr,ijitter,iampJitter,:,r) = taskAmp*(1+ampJitter(iampJitter)*randn(1))*runTC(1:upsampleFactor:end);%downsample
                rwdSignal(isnr,ijitter,iampJitter,:,r) = runTC(1:upsampleFactor:end);%downsample
       
            end
            n = (taskAmp/snr(isnr))*randn(size(rwdSignal(isnr,ijitter,iampJitter,:,:)));
            
            rwdNoise(isnr,ijitter,iampJitter,:,:) = abs(ifft(repmat(oneOverF,1,1,1,1,runsPerRwd).*fft(n)));
        end
    end
end

%%
rwdTC = rwdSignal + rwdNoise;
rwdTC = zscore(rwdTC,0,4);
rwdTC = rwdTC(:,:,:,trialLength+1:end,:);%junk first cycle
for isnr=1:length(snr)
    for ijitter = 1:length(jitter)
        for iampJitter=1:length(ampJitter)
            for r=1:runsPerRwd
                rwdTC(isnr,ijitter,iampJitter,:,r) = filterTC(squeeze(rwdTC(isnr,ijitter,iampJitter,:,r)));
            end
        end
    end
end



rwdTrials = reshape(rwdTC,length(snr),length(jitter),length(ampJitter),trialLength,[]);
meanTrial = mean(rwdTrials,5);
meanRun = squeeze(mean(rwdTC,5));
for isnr=1:length(snr)
    for ijitter = 1:length(jitter)
        for iampJitter=1:length(ampJitter)
            f = squeeze(fft(rwdTrials(isnr,ijitter,iampJitter,:,:)));
            fftTrialAmp(isnr,ijitter,iampJitter,:) = abs(f(2,:));
            fftTrialPh(isnr,ijitter,iampJitter,:) = angle(f(2,:));
            fftRun(isnr,ijitter,iampJitter,:) = abs(fft(meanRun(isnr,ijitter,iampJitter,:)));
        end
    end
end

nframes = length(fftRun);
%from plotMeanFourierAmp.m
fftRun = fftRun / (nframes/2);
frequencies = [0:nframes-1]/(nframes*TR);
fftRun = fftRun(:,:,:,1:floor(nframes/2));
frequencies = frequencies(1:floor(nframes/2));



%%
trialsAmp = squeeze(std(rwdTrials,0,4));
ampMeanTrial = std(meanTrial,0,4);

meanTrialsAmp = mean(trialsAmp,4);
stdTrialsAmp = std(trialsAmp,0,4);

meanTrialsFFTamp = mean(fftTrialAmp,4);

phVar = circ_std(fftTrialPh,[],[],4);
ampVar = std(fftTrialAmp,0,4);
totalVar = mean(std(rwdTrials,0,5),4);



%% PLOT VARIABILITY AS FUNCTION OF SNR AND TEMPORAL JITTER
i=i+1; figure(i) ;clf
rows=3;
cols=7;
xgv = 1:size(ampMeanTrial,1);
ygv = 1:size(ampMeanTrial,2);
zgv = 1:size(ampMeanTrial,3);
[X,Y,Z] = meshgrid(xgv,ygv,zgv);
for r=1:rows
    xsize = npoints;
    ysize = npoints;
    switch r
        case 1
            ind =  Z==1;
            xlabelstring = 'temporal jitter (radians)';
            ylabelstring = 'SNR';
        case 2
            ind =  X==1;
            xlabelstring = 'amplitude jitter';
            ylabelstring = 'SNR';
        case 3
            ind =  Y==1;
            xlabelstring = 'amplitude jitter';
            ylabelstring = 'temporal jitter (radians)';
    end
    c=1;
    subplot(rows,cols,(r-1)*cols+c)
    temp = ampMeanTrial(ind);
    temp = reshape(temp,xsize,ysize);
    imagesc(temp)
    title('amplitude of mean trial');
    ylabel(ylabelstring);
    xlabel(xlabelstring);
    
    c=c+1;
    subplot(rows,cols,(r-1)*cols+c)
    temp = meanTrialsAmp(ind);
    temp = reshape(temp,xsize,ysize);
    imagesc(temp)
    title('mean trials amplitude');
    ylabel(ylabelstring);
    xlabel(xlabelstring);
    
        c=c+1;
    subplot(rows,cols,(r-1)*cols+c)
    temp = meanTrialsFFTamp(ind);
    temp = reshape(temp,xsize,ysize);
    imagesc(temp)
    title('mean trials FFT amplitude');
    ylabel(ylabelstring);
    xlabel(xlabelstring);
    
    
    c=c+1;
    subplot(rows,cols,(r-1)*cols+c)
    temp = stdTrialsAmp(ind);
    temp = reshape(temp,xsize,ysize);
    imagesc(temp)
    title('amplitude variability');
    ylabel(ylabelstring);
    xlabel(xlabelstring);
    
    c=c+1;
    subplot(rows,cols,(r-1)*cols+c)
    temp = ampVar(ind);
    temp = reshape(temp,xsize,ysize);
    imagesc(temp)
    title('FFT amp variability');
    ylabel(ylabelstring);
    xlabel(xlabelstring);
    

    
    
    
    c=c+1;
    subplot(rows,cols,(r-1)*cols+c)
    temp = phVar(ind);
    temp = reshape(temp,xsize,ysize);
    imagesc(temp)
    title('FFT phase variability');
    ylabel(ylabelstring);
    xlabel(xlabelstring);
    
    c=c+1;
    subplot(rows,cols,(r-1)*cols+c)
    temp = totalVar(ind);
    temp = reshape(temp,xsize,ysize);
    imagesc(temp)
    title('mean timepoint variability');
    ylabel(ylabelstring);
    xlabel(xlabelstring);
end


%% 
linewidth = 0.1;
markersize = 10;
i=i+1; figure(i) ;clf
rows=3; cols=3;

%SNR
subplot(rows,cols,1)
ijitter=1;
iampJitter=1;
for isnr=1:length(snr)    
    plot(squeeze(meanTrial(isnr,ijitter,iampJitter,:)),'color',plotColorsSnr{isnr},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,2)
for isnr=1:length(snr)
    plot(frequencies,squeeze(fftRun(isnr,ijitter,iampJitter,:)),'.-','color',plotColorsSnr{isnr},'linewidth',linewidth,'markersize',markersize);
%     plot(frequencies,log(squeeze(fftRun(isnr,ijitter,:))),'.-','color',plotColorsSnr{isnr},'linewidth',linewidth,'markersize',markersize);
    hold on
end
% ylabel('log(FFT amplitude)');
subplot(rows,cols,3)
for isnr=1:length(snr)
    plot(squeeze(std(rwdTrials(isnr,ijitter,iampJitter,:,:),0,5)),'color',plotColorsSnr{isnr},'linewidth',linewidth);
    hold on
end

%temporal variability
subplot(rows,cols,cols+1)
isnr=1;
for ijitter=1:length(jitter)    
    plot(squeeze(meanTrial(isnr,ijitter,iampJitter,:)),'color',plotColorsJitter{ijitter},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,cols+2)
for ijitter=1:length(jitter)
    plot(frequencies,squeeze(fftRun(isnr,ijitter,iampJitter,:)),'.-','color',plotColorsJitter{ijitter},'linewidth',linewidth,'markersize',markersize);
%     plot(frequencies,log(squeeze(fftRun(isnr,ijitter,:))),'.-','color',plotColorsJitter{ijitter},'linewidth',linewidth,'markersize',markersize);
    hold on
end
% ylabel('log(FFT amplitude)');
subplot(rows,cols,cols+3)
for ijitter=1:length(jitter)
    plot(squeeze(std(rwdTrials(isnr,ijitter,iampJitter,:,:),0,5)),'color',plotColorsJitter{ijitter},'linewidth',linewidth);
    hold on
end

%amplitude variability
subplot(rows,cols,2*cols+1)
isnr=1;
ijitter=1;
for iampJitter=1:length(ampJitter)    
    plot(squeeze(meanTrial(isnr,ijitter,iampJitter,:)),'color',plotColorsAmpJitter{iampJitter},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,2*cols+2)
for iampJitter=1:length(ampJitter)
    plot(frequencies,squeeze(fftRun(isnr,ijitter,iampJitter,:)),'.-','color',plotColorsAmpJitter{iampJitter},'linewidth',linewidth,'markersize',markersize);
%     plot(frequencies,log(squeeze(fftRun(isnr,ijitter,:))),'.-','color',plotColorsJitter{ijitter},'linewidth',linewidth,'markersize',markersize);
    hold on
end
% ylabel('log(FFT amplitude)');
subplot(rows,cols,2*cols+3)
for iampJitter=1:length(ampJitter)
    plot(squeeze(std(rwdTrials(isnr,ijitter,iampJitter,:,:),0,5)),'color',plotColorsAmpJitter{iampJitter},'linewidth',linewidth);
    hold on
end

keyboard
%%
i=i+1; figure(i) ;clf
rows=3; cols=2;
linewidth = 1;
markersize = 10;
%SNR
subplot(rows,cols,1)
ijitter=1; iampJitter=1;
for isnr=1:length(snr)    
    plot(TR*(0:trialLength-1),squeeze(meanTrial(isnr,ijitter,iampJitter,:)),'color',plotColorsSnr{isnr},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,2)
for isnr=1:length(snr)
    plot(TR*(0:trialLength-1),squeeze(std(rwdTrials(isnr,ijitter,iampJitter,:,:),0,5)),'color',plotColorsSnr{isnr},'linewidth',linewidth);
    hold on
end

%temporal variability
subplot(rows,cols,cols+1)
isnr=1;
for ijitter=1:length(jitter)    
    plot(TR*(0:trialLength-1),squeeze(meanTrial(isnr,ijitter,iampJitter,:)),'color',plotColorsJitter{ijitter},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,cols+2)
for ijitter=1:length(jitter)
    plot(TR*(0:trialLength-1),squeeze(std(rwdTrials(isnr,ijitter,iampJitter,:,:),0,5)),'color',plotColorsJitter{ijitter},'linewidth',linewidth);
    hold on
end

%amplitude variability
subplot(rows,cols,2*cols+1)
isnr=1;
ijitter=1;
for iampJitter=1:length(ampJitter)    
    plot(TR*(0:trialLength-1),squeeze(meanTrial(isnr,ijitter,iampJitter,:)),'color',plotColorsAmpJitter{iampJitter},'linewidth',linewidth);
    hold on
end
subplot(rows,cols,2*cols+2)
for iampJitter=1:length(ampJitter)
    plot(TR*(0:trialLength-1),squeeze(std(rwdTrials(isnr,ijitter,iampJitter,:,:),0,5)),'color',plotColorsAmpJitter{iampJitter},'linewidth',linewidth);
    hold on
end

for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    xlabel('time (sec)');
    if mod(isubplot,2)==1
        ylabel('BOLD signal (std image intensity)');
    else
        ylabel('signal variability (std)');
    end
    drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -8/64, 'xAxisMargin', 4/64, 'yAxisMargin', 0/64,'xAxisMinMaxSetByTicks',0,...
        'labelFontSize',7);
    axis square
end
set(gcf,'position',[10 10 15 24]);
print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig8.pdf']);


%%
% i=i+1; figure(i) ;clf
% rows=1;
% cols=npoints;
% for i=1:npoints
%     subplot(rows,cols,i)
%     histogram(taskAmp*(1+atan(ampJitter(i)*randn(1000,1))))
%     hold all
% end
%%
% i=i+1; figure(i) ;clf
% rows=1;
% cols=2;
% subplot(rows,cols,1)
% plot(squeeze(rwdTC(1,1,1,:,1)))
% subplot(rows,cols,2)
% plot(squeeze(rwdTC(1,end,1,:,1)))

%%
function filteredTC = filterTC(tc)
filt = ones(1,length(tc));
filtNonzeros = [0,0,0.117503097415405,0.588887709492812,0.904365555167461,0.988891003461758,0.999355620179431];
filt(2:8) = filtNonzeros;
filt(end-6:end) = filtNonzeros;
filteredTC = real(ifft(fft(tc) .* filt' ))';

end