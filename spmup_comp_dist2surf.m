function davg = spmup_comp_dist2surf(anat)

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
% Remi Gau - University of Birmingham


fprintf('Computing average distance to brain surface.\n')

[path, file] = spm_fileparts(anat);

try

surface_file = spm_select('FPList', path, ['^c1' file '.*\.surf\.gii$']);

if isempty(surface_file)
    
    fprintf(' No brain surface was found. Trying to create one.\n')
    % create a surface of the brain using the TPM if we don't have one
    files = spm_select('FPList', path, ['^c[12]' file '.*\.nii$']);
    if isempty(files)
        warning(' Could not find TPMs resulting from brain segmentation.\n')
    end
    spm_surf(files, 2);
    surface_file = spm_select('FPList', path, ['^c1' file '.*\.surf\.gii$']);
    
end

% compute the average surface to the brain
FV = gifti(surface_file);
center = FV.vertices(FV.faces(:, :), :);
center = reshape(center, [size(FV.faces,1) 3 3]);
center = squeeze(mean(center,2));
ori_dist = sqrt(sum((center.*-1).^2,2))';
davg = mean(ori_dist);

catch

    fprintf(' Could not compute the average distance to the brain surface.\n')
    fprintf(' Using the default value instead.\n')
    davg = 50; % we give it a default value if anything fails

end

fprintf(' Average distance to the cortex surface: %f mm \n', davg)


end

