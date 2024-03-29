function options = spmup_getoptions(BIDS_dir)

% routine to return all the spmup options with defaults set
%
% FORMAT options = spmup_getoptions(BIDS_dir)
% INPUT BIDS_dir is a folder in which analyses should be carried on
% OUTPUT options is a structure with spmup optional default parameters set
%               'outdir' where to write the data (default is [BIDS_dir filesep 'derivatives'])
%               'keep_data' is 'off' (default) or 'on' to keep all steps 
%                           off means only the last processed data are available
%               'overwrite_data' 'on' (default) turning it 'off' is useful to restart
%               	             some processing while kepping previous steps
%               'Ncores' [] default or an integer to be used by the matlab parallel processing toolbox
%                        if left empty it will set to use Ncores-1;
%               'QC' 'on' (default) or 'off' performs quality control (also augment motion regressor file) 
%               'subjects' a cell array of subject subset from BIDS dataset (if empty, default, do them all)
%               'anat' a cell array of file pattern to use e.g. {'acq_MPRAGE_T1w','acq_SPACE_T2w'} 
%               'reorient' is 'on' (default) to automatically set (0,0,0)
%               'removeNvol' set to '0' (default) is the number of initial volumes to remove
%               'ses'  [] (default) to specificy which session to analyze.
%               'task' [] (default) to specificy which bold task to analyze.
%               'acq'  [] (default) to specificy which bold acquisition to analyze.
%               'rec'  [] (default) to specificy which bold reconstruction to analyze.
%               'despike' can be 'before', 'after' or 'off' (default) runs median despiking after or before realignment 
%               'fmap' is 'on' (default) to indicate to use field maps (can be set to 'off')
%                      this can also be set to one of BIDS valid field when
%                      multiple types maps are available e.g. 'phasediff' or 'epi'
%               'VDM' cellarray of file on drive (to apply field map distorsion correction - see spmup_compute_vdm.m)
%                     this option takes precedence over 'fmap' is not empty (default)
%               'multiecho' 'no' (default) or 'yes' denotes whether bold are multi-echo
%               'norm' 'T1norm' (default) or 'EPInorm' choice of the type of template for normalization
%               'norm_res' can be any voxel size, 'EPI' or 'EPI-iso' (default - round to nearest isotropic voxel size)
%               'norm_inv_save', is you want to save the 'big' inverse transform file 'off' by default
%               'skernel' is the smoothing kernel as how many times the voxel size
%                         [] is the default because 'derivatives' default is 1 (and hrf re-estimation needs re-smoothing) 
%                         using [] implies smooth by norm_res, GLM, correct, and re-smooth 3*norm_res (task fmri)
%                         this can also be set to 'off' to skip smoothing altogether
%               'motionexp' 'off' (default) or 'on' compute voltera motion expansion
%               'scrubbing' 'on' (default) or 'off' find outliers in motion derivatives and in globals
%               'nuisance' 'WMCSF' (default) or 'compcor' for cleaning resting state data
%               'GLM', 'on' (default) or 'off' to compute the simple GLM (see spm_BIDS_1rstlevel.m)
%               'highpass' [] by default means 128 sec if task or 1000 sec for rest
%                'conditions', 'trial_type' (defaut) is the column name from tsv file to build default GLM cndition regressors from
%               'derivatives' to include '1'st (default), '2'nd or '0' derivatives in the model
%                             --> when the 1st derivative is included, the hrf amplitude estimates are
%                                 automatcally 'boosted' and time delay also computed
%               'carpet_plot' is 'on' (default) or 'off' to display GLM residuals voxels in 2D 
%                            (like for resting state - if denoised, and regressed propertly, 
%                             residuals should be clean)
%
% see also see spmup_compute_vdm.m and spm_BIDS_1rstlevel.m 
%
% Cyril Pernet July 2022
% --------------------------
%  Copyright (C) SPMUP Team 

if strcmpi(BIDS_dir(end),filesep)
    BIDS_dir = BIDS_dir(1:end-1);
end

options = struct(...
    'outdir', [BIDS_dir filesep 'derivatives'], ...
    'keep_data', 'off',  ...
    'overwrite_data', 'on', ...
    'Ncores',[],...
    'QC', 'on', ...
    'subjects',[], ...
    'anat',[], ...
    'reorient','on', ...
    'ses',[], ...
    'task',[], ...
    'acq',[], ...
    'rec',[], ...
    'removeNvol', 0, ...
    'despike', [], ...
    'fmap','on',...
    'VDM', [], ...
    'multiecho', 'no', ...
    'norm', 'T1norm', ...
    'norm_res', 'EPI-iso', ...
    'norm_inv_save','off',...
    'skernel', [], ...
    'motionexp', 'off', ...
    'scrubbing', 'on', ...
    'nuisance','WMCSF',...
    'GLM', 'on', ...
    'highpass', [], ...
    'conditions', 'trial_type', ...
    'derivatives', 1, ...
    'carpet_plot', 'on');
