
curFolder = pwd;
% dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects_concatenated/';
% subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', ...
%     '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', ...
%     '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};

dataFolder = '/Volumes/MH02086153MACDT-Drobo/allMinSubjects/';
subFolders = {'000520180116', '0008i20180213', '0016i20180207', '002220171212', '003220180105', '0034i20180209', '003520180328', '004020180328','004120180320', '0042i20180412', '0045i20180309', '0046i20180409', '0049i20180404', '005220180621'};
numSubs=length(subFolders);
cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end
trialLength=10;
clear concatInfo  runRwd

rows=2;
cols=ceil(numSubs/rows);
for iSub = 1:numSubs
    
    cd(subFolders{iSub});
    v=newView;
    
    % Find the correct Concatenation group
    concatGroupNum = viewGet(v,'groupNum','Concatenation'); %concatenation of multiple sessions
    if isempty(concatGroupNum)
        iSub
%         concatGroupNum = viewGet(v,'groupNum','Concatenation');
    end
    
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', concatGroupNum);
    nScans = viewGet(v, 'nscans');
    for iScan = 1:nScans%2 concatenations, 1 for each reward type
        s = viewGet(v, 'stimfile', iScan);
        %             if ~isfield(s{1}, 'stimulus')
        %                 s=s{1};
        %             end

        concatInfo{iSub,iScan} = viewGet(v, 'concatInfo', iScan);
        for r=1:length(s)
            runRwd{iSub,iScan}(r) = s{r}.stimulus.rewardVal;
        end
        
    end
    
    deleteView(v);
    cd ..
end
set(gcf,'position',[100 100 1000 500]);
%%


runRwd

