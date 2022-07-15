function davg = spmup_comp_dist2surf(anat,c1,c2)

% Compute the average surface to the brain. Will return a default value of
% 50 mm if this fails.
%
% Adapted from motion finger print functions and script (mw_mfp.m) from Marko Wilke
% https://www.medizin.uni-tuebingen.de/kinder/en/research/neuroimaging/software/
% see http://www.dx.doi.org/10.1016/j.neuroimage.2011.10.043
% and http://www.dx.doi.org/10.1371/journal.pone.0106498
%
% FORMAT: davg = spmup_comp_dist2surf(anat)
%
% INPUT: anat is the fullpath to the T1w image used to compute the distance
%           to the brain surface. It must have been segmented by SPM and the
%           function will look for the c1* and c2* of the grey and white
%           matter TPMs.
%
% OUTPUT: davg is the average distance to the brain surface
%
% Remi Gau & Cyril Pernet
% --------------------------
%  Copyright (C) SPMUP Team 

davg = [];

%% deal with inputs

[path, file] = spm_fileparts(anat);
if ~exist(anat,'file')
    error(' anat file %s doesn''t exist\n',file);
end

if nargin == 1
    surface_file = [];
    surface_file = spm_select('FPList', path, ['^c1' file '.*\.surf\.gii$']);
else
    spath        = spm_fileparts(anat);
    files(1,:)   = c1;
    files(2,:)   = c2; 
    spm_surf(files, 2);
    surface_file = spm_select('FPList', spath, '.*\.surf\.gii$');
end

if isempty(surface_file)
    fprintf(' No brain surface was found. Trying to create one.\n')
    files = spm_select('FPList', path, ['^c[12]' file '.*\.nii$']);
    if isempty(files)
        warning(' Could not find c1 c2 files resulting from brain segmentation.\n')
        warning(' Returning the default value 50 instead.\n')
        davg = 50; % we give it a default value if anything fails
    else
        spm_surf(files, 2);
        surface_file = spm_select('FPList', spath, '.*\.surf\.gii$');
    end    
end

%% compute the average surface to the brain
if isempty(davg)
    if isempty(surface_file)
        warning(' Could not find surface gii file from runing spm_surf.\n')
        warning(' Returning the default value 50 instead.\n')
        davg = 50; % we give it a default value if anything fails
    else
        FV       = gifti(surface_file);
        center   = FV.vertices(FV.faces(:, :), :);
        center   = reshape(center, [size(FV.faces,1) 3 3]);
        center   = squeeze(mean(center,2));
        ori_dist = sqrt(sum((center.*-1).^2,2));
        davg     = mean(ori_dist);
        fprintf(' Average distance to the cortex surface: %f mm \n', davg)
    end
end
