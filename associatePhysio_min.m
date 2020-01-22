%stimfiles must be already associated
clear d d2 d3 
% newPwd = fullfile('/Users/rothzn/min_all_subjects',curSubname);
v=newView;

v = viewSet(v,'groupname','Raw');
nScans = viewGet(v, 'nscans');

[FILEPATH,NAME,EXT] = fileparts(pwd);
pwdsplit = split(FILEPATH,'/');
curSubname = NAME;

% origFullfile = fullfile(newPwd,'origFilenames.mat');

% %find realtime directory
% d=dir(fullfile(newPwd, '*.tar'));
% d2 = dir(fullfile(newPwd,[d.name(1:12) '*']));%beginning of the name of the zipped data file
% d3 = dir(fullfile(newPwd,d2(1).name));%should include only one folder
% realtimePath = fullfile(newPwd, d2(1).name, d3(end).name, 'realtime/');

%find realtime directory
realtimePath = fullfile(pwd, 'Etc');

%go through scans
for iscan = 1:nScans
    s = viewGet(v, 'stimfile', iscan);
    filename = viewGet(v,'tseriesfile', iscan);%should be corrected filename
    
    k = strfind(filename,'mr_00');
    if length(k)==1
        numrun = str2num(filename(k+5:k+6));
        numrun = numrun+1;
        numrunstr = num2str(numrun,'%04u');
    else%uncorrected filename?
        keyboard
    end
    
    
%     d = dir([realtimePath '*scan_' numrunstr '*']);%d(1) is ecg, d(2) is resp
    d = dir(fullfile(realtimePath, ['*scan_' numrunstr '*']));%d(1) is ecg, d(2) is resp
    s{1}.myscreen.originalmlrdir = pwd;
    s{1}.myscreen.physiodir = realtimePath;
    s{1}.myscreen.ecgfilename = d(1).name;
    s{1}.myscreen.respirationfilename = d(2).name;
    stimfilename = viewGet(v,'stimfilename',iscan);
    [FILEPATH,NAME,EXT] = fileparts(stimfilename{1});
    %     save(stimfilename{1},s{1}.myscreen, s{1}.task, s{1}.
    temp = s{1};
    %         save(['~/temp/' [NAME EXT]],'-struct','temp');
    save(stimfilename{1},'-struct','temp');
end
 
deleteView(v);