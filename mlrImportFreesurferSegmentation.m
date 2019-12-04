function[] = mlrImportFreesurferSegmentation()
% mlrImportFSLabels
%
% Converts freesurfer cerebellum segmentation into mlr coordinates
%
% znr 11/2019

aseg = MRIread('mri/aseg.mgz');
%FreeSurferColorLut.txt from /Applications/freeesurfer/
% labels: 8 - left cerebellum cortex, 47 - right cerebellum cortex
labels = [8 47];
labelNames = {'leftCerebellarCortex','rightCerebellarCortex'};
% grab the xform from the base anatomy
anatFileName = dir('surfRelax/*mprage_pp.hdr');
h = cbiReadNiftiHeader(fullfile('surfRelax', anatFileName.name));

for iLabel = 1:length(labels)
    roi.voxelSize = h.pixdim(2:4);
    % set the xform and xformCode to whatever is in the curent scan
    roi.xform = h.qform44;
    roi.sformCode = 0;
    volumeCropSize = h.dim(2:4)';
    
    roi.name = labelNames{iLabel};
    
    volSize = size(aseg.vol);
    ind = find(aseg.vol==labels(iLabel));
    clear coords
    [coords(:,1), coords(:,2), coords(:,3)] = ind2sub(volSize,ind);
%     convertedCoords = Coords(:,2:4)+volumeCropSize./2 + 1;



% volumeCropSize = volSize;
% coords = volSize - coords;
% coords(:,3) = volSize(3) - coords(:,3);

coords(:,1) = volSize(1) - coords(:,1);
coords(:,2) = volSize(2) - coords(:,2);
% coords(:,2) = coords(:,2) - 40;

%     convertedCoords = [coords(:,2)*volumeCropSize(1)/volSize(1) ...
%     coords(:,3)*volumeCropSize(2)/volSize(2) ...
%     coords(:,1)*volumeCropSize(3)/volSize(3)];
    
%     convertedCoords = coords;
    convertedCoords = [coords(:,2) coords(:,3) coords(:,1)] - (volSize-volumeCropSize)./2;
    convertedCoords = convertedCoords + [0 0 0];
    
    convertedCoords(:,4) = 1;

    roi.coords = convertedCoords';
    % make it into an roi
    [tf roi] = isroi(roi);
    %save roi as a cell array
    
    
    % path to file
    roiName = fixBadChars(stripext(labelNames{iLabel}),[], {'.','_'});
    eval([roiName,'=roi;']);
    filename = roiName;
    pathStr = fullfile('surfRelax',filename);
    saveString = ['save(pathStr,','''',roiName,'''',');'];
    eval(saveString);
    
%     keyboard
end
keyboard
