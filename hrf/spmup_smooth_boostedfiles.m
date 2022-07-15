function out = spmup_smooth_boostedfiles(P,K)

% smooth the data in P with a Gaussian kernel of FHWM [x y z]
% and masks out using the mask located in folder above your boosted files
%
% FORMAT: out = spmup_smooth_boostedfiles(P,[x y z]);
%
% INPUT: P names of files to be smmothed
%        x y z the FWHM of the Gaussian kernel
%
% OUTPUT: out the names of the smoothed files
%
% see also spm_smooth
% Cyril Pernet March 2016
% ----------------------------
% Copyright (c) SPMU+ toolbox

for file = size(P,1):-1:1
    % get the mask
    [sp,sf,se] = fileparts(P(file,:));
    V(1)       = spm_vol([fileparts(sp) filesep 'mask.nii']);
    mask       = spm_read_vols(V(1)) == 0; % make sure it's binary/logical
    % smooth
    spm_smooth(P(file,:), [sp filesep 's' sf se],K,0)
    % mask out
    V(2)       = spm_vol([sp filesep 's' sf se]);
    data       = spm_read_vols(V(2));
    data(mask) = NaN;
    V          = spm_write_vol(V(2),data);
    out{file} = V.fname;
end
