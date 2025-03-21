function [rho, pval, regionList, xRegional, yRegional] = regionalCorr(atlasFile, xFiles, yFile, corrMethod)
% regionalCorr computes correlations between regional cell proportion values
% (extracted from a list of NIfTI files, one per cell type) and a neuroimaging phenotype.
%
% INPUTS:
%   atlasFile - string; filename of the atlas NIfTI file (defines regions)
%   xFiles    - cell array of strings; filenames of NIfTI files for cell proportions.
%               Each file corresponds to one cell type (e.g., {'ast.nii','end.nii', 'mic.nii', 'neu.nii', 'oli.nii', 'opc.nii'})
%   yFile     - string; filename of the neuroimaging phenotype NIfTI file of interest (e.g., atrophy map)
%   corrMethod- (optional) string; type of correlation to compute, 'Spearman' (default) or 'Pearson'
%
% OUTPUTS:
%   rho       - vector; correlation coefficients (one per cell type)
%   pval      - vector; p-values associated with the correlations
%   regionList- vector; unique region identifiers as defined in the atlas (after removal of regions with 0 or NaN)
%   xRegional - matrix; extracted average cell proportion values per region for each cell type
%               (dimensions: number of regions x number of cell types)
%   yRegional - vector; extracted average neuroimaging phenotype value per region
%
% EXAMPLE USAGE:
%   atlasFile = 'DTKatlas.nii';
%   xFiles = {'ast.nii', 'end.nii', 'mic.nii', 'neu.nii', 'oli.nii', 'opc.nii'};
%   yFile = 'LOAD_atrophy.nii';
%   [rho, pval, regionList, xRegional, yRegional] = cellNeuroCorrFiles(atlasFile, xFiles, yFile, 'Pearson');
%
% NOTE:
%   This function uses SPM (e.g., spm_vol) to obtain affine transformation matrices.
%   Ensure that SPM is added to your MATLAB path.

%% Set default correlation method if not provided
if nargin < 4
    corrMethod = 'Spearman';
end

%% Step 1: Read the atlas and extract its region information
atlas = niftiread(atlasFile);
ind = find(atlas);              % indices of nonzero (i.e., defined) voxels
regionLabels = atlas(ind);      % region labels at these voxels

% Convert linear indices to voxel coordinates in the atlas
[row, col, zed] = ind2sub(size(atlas), ind);
voxelMask = [row, col, zed];

% Get unique region identifiers from the atlas
regionList = unique(regionLabels);

% Load the atlas volume using SPM to obtain the affine transformation matrix
atlasVol = spm_vol(atlasFile);
mniMask = atlasVol.mat * [voxelMask, ones(size(voxelMask,1),1)]';
mniMask = mniMask';
mniMask(:,4) = [];  % remove homogeneous coordinate
mniMask = round(mniMask);

%% Step 2: Extract regional values from cell proportion (X) files
nCellTypes = length(xFiles);
nRegions   = length(regionList);
xRegional = nan(nRegions, nCellTypes);

for i = 1:nCellTypes
    xRegional(:, i) = extractRegionalFromFile(xFiles{i}, mniMask, regionLabels, regionList);
end

%% Step 3: Extract regional values from the neuroimaging phenotype (Y)
yRegional = extractRegionalFromFile(yFile, mniMask, regionLabels, regionList);

%% Step 4: Remove regions (rows) with any 0 or NaN in X or Y, and report the row numbers
rowsToRemove = find( any(xRegional == 0 | isnan(xRegional), 2) | (yRegional == 0 | isnan(yRegional)) );
if ~isempty(rowsToRemove)
    fprintf('Removing %d region(s) due to a 0 or NaN value in X or Y. Row number(s): %s\n', ...
            numel(rowsToRemove), mat2str(rowsToRemove));
    regionList(rowsToRemove) = [];
    xRegional(rowsToRemove, :) = [];
    yRegional(rowsToRemove) = [];
end

%% Step 5: Compute correlation for each cell type across regions
rho  = nan(nCellTypes, 1);
pval = nan(nCellTypes, 1);

for i = 1:nCellTypes
    [r, p] = corr(xRegional(:, i), yRegional, 'Type', corrMethod, 'Rows', 'complete');
    rho(i)  = r;
    pval(i) = p;
end

end

%% Local helper function to extract average regional values from a NIfTI file
function regVals = extractRegionalFromFile(niftiFile, mniMask, regionLabels, regionList)
    % Read the NIfTI file
    map = niftiread(niftiFile);
    
    % Find indices of nonzero voxels in the file
    ind2 = find(map);
    [row2, col2, zed2] = ind2sub(size(map), ind2);
    voxels = [row2, col2, zed2];
    
    % Load the volume using SPM to obtain its affine transformation matrix
    vol = spm_vol(niftiFile);
    if iscell(vol)
        vol = vol{1};
    end
    mniMap = vol.mat * [voxels, ones(size(voxels,1),1)]';
    mniMap = mniMap';
    mniMap(:,4) = [];  % remove homogeneous coordinate
    mniMap = round(mniMap);
    
    nRegions = length(regionList);
    regVals = nan(nRegions, 1);
    
    % Loop through each region defined in the atlas
    for j = 1:nRegions
        % Find indices in the atlas corresponding to the current region
        regionInds = find(regionLabels == regionList(j));
        if isempty(regionInds)
            regVals(j) = NaN;
            continue;
        end
        
        % Get the MNI coordinates for voxels in the current region (from the atlas)
        regionCoords = mniMask(regionInds, :);
        
        % Find matching voxels in the current NIfTI file based on MNI coordinates
        [~, loc] = ismember(regionCoords, mniMap, 'rows');
        loc = nonzeros(loc);  % remove zeros (no match)
        
        if isempty(loc)
            regVals(j) = NaN;
        else
            regionValues = zeros(length(loc), 1);
            for k = 1:length(loc)
                coord = voxels(loc(k), :);
                regionValues(k) = map(coord(1), coord(2), coord(3));
            end
            % Average the voxel values in this region
            regVals(j) = mean(regionValues);
        end
    end
end
