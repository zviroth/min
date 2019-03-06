
dataFolder = '/Volumes/MH02086153MACDT-Drobo/decodingAnalysis/rwd/s004520180309/Etc/';
dataFolder = '/Volumes/MH02086153MACDT-Drobo/decodingAnalysis/rwd/s005220180621/Etc/';
cd(dataFolder);
clear ecg resp

%% ECG
d=dir('ECG*.*');

figure(1)
clf
rows = ceil(length(d)/3);
cols = ceil(length(d)/rows);
select=0.03;
display=0;
for i=1:length(d)
    filename = d(i).name;
    ecg{i}=load(filename);
    subplot(rows,cols,i)
    plot(ecg{i})
    hold on
    [ecgPeaks{i},criterion] = pickpeaks(ecg{i},select,display);
    
    ecgPeaksDiff{i} = diff(ecgPeaks{i});
    ecgPeaksAmp{i} = ecg{i}(ecgPeaks{i});
    plot(ecgPeaks{i},ecgPeaksAmp{i},'ro','MarkerSize',3,'LineWidth',2)

end

figure(2)
clf
for i=1:length(d)
    
    subplot(rows,cols,i)
    plot(ecgPeaksDiff{i})
    hold all
%     plot(ecgPeaksAmp{i})
end



%% RESPIRATION

d=dir('Resp*.*');

figure(3)
clf
rows = ceil(length(d)/3);
cols = ceil(length(d)/rows);
select=0.1;
display=0;
for i=1:length(d)
    filename = d(i).name;
    resp{i}=load(filename);
    subplot(rows,cols,i)
    plot(resp{i})
    hold on
    [respPeaks{i},criterion] = pickpeaks(resp{i},select,display);
    plot(respPeaks{i},resp{i}(respPeaks{i}),'ro','MarkerSize',3,'LineWidth',2)
    respPeaksDiff{i} = diff(respPeaks{i});
    respPeaksAmp{i} = ecg{i}(respPeaks{i});
%     plot(respPeaks{i},respPeaksAmp{i},'ro','MarkerSize',3,'LineWidth',2)
end

figure(4)
clf
for i=1:length(d)
    
    subplot(rows,cols,i)
    plot(respPeaksDiff{i})
    hold all
%     plot(respPeaksAmp{i})
end
