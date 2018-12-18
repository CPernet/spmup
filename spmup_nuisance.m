function nuisances = spmup_nuisance(fmridata,whitematter,csf)

% routine to compute nuisance regressors from White matter and CSF masks
% masks are thresholded at 70%, the time series extracted, detrended (1st
% order polynomial) and averaged
%
% FORMAT: nuisances = spmup_nuisance(fmridata,whitematter,csf)
%
% IMPUT fmridata is the time series of data
%       whitematter is the WM mask
%       csf is the CSF mask
%
% OUTPUT nuisances is a structure with fields .WM and .CSF
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

threshold = 0.7;

if nargin == 0
    [fmridata,sts]= spm_select(Inf,'image','Select fMRI data');
    if sts == 0; disp('selection aborded'); return; end
    if size(fmridata,1) == 1 && strcmp(fmridata(length(fmridata)-1:end),',1')
        fmridata = fmridata(1:length(fmridata)-2); % in case picked 4D put left ,1
    end
    [whitematter,sts]= spm_select(1,'image','Select white matter tissue image');
    if sts == 0; disp('selection aborded'); return; end
    [csf,sts]= spm_select(1,'image','Select csf tissue image');
    if sts == 0; disp('selection aborded'); return; end
end

% check fmridata is 4D
if iscell(fmridata)
    for v=1:size(fmridata,1)
        Vfmri(v) =spm_vol(fmridata{v});
    end
else
    Vfmri = spm_vol(fmridata);
end

if sum(size(Vfmri)) == 2
    error('fMRI data must be time series')
end
    
Vwhitematter = spm_vol(whitematter);
if any(Vwhitematter.dim ~= Vfmri(1).dim)
    error('Dimension issue between data and white matter')
end
whitematter = spm_read_vols(Vwhitematter);

Vcsf = spm_vol(csf); 
if any(Vcsf.dim ~= Vfmri(1).dim)
    error('Dimension issue between data and CSF')
end
csf = spm_read_vols(Vcsf);

% check if the masks are binary with a threshold 
whitematter = whitematter > threshold ;
csf = csf > threshold ;

% noise regressor comes from the white matter and csf images
index = find(whitematter);
[x,y,z]=ind2sub(size(whitematter),index);
nuisances.WM = nanmean(spm_detrend(spm_get_data(Vfmri,[x y z]'),1),2);

index = find(csf);
[x,y,z]=ind2sub(size(csf),index);
nuisances.CSF = nanmean(spm_detrend(spm_get_data(Vfmri,[x y z]'),1),2);
