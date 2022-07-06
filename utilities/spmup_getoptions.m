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
%               'QC' 'on' (default) or 'off' performs quality control (also augment motion regressor file) 
%               'anat' a cell array of file pattern to use e.g. {'_T1w','_T2w'} 
%               'removeNvol' set to '0' (default) is the number of initial volumes to remove
%               'VDM' cellarray (either a file on drive or a memory mapped variable)
%               'task' [] (default) to specificy which bold task to analyze.
%               'rec' [] (default) to specificy which bold reconstruction to analyze. 
%               'acq' [] (default) to specificy which bold acquisition to analyze. 
%               'despike' 'on' (default) or 'off' runs median despiking
%               'motionexp' 'off' (default) or 'on' compute voltera motion expansion
%               'scrubbing' 'on' (default) or 'off' find outliers in motion derivatives and in globals
%                'norm' 'T1norm' (default) or 'EPInorm' choice of the type of template for normalization
%               'ignore_fieldmaps' 'on' or 'off' (default) to ignore distorsion correction for T1norm
%               'skernel' 3 by default is the smoothing kernel as how many times the voxel size
%
% Cyril Pernet July 2022
% --------------------------
%  Copyright (C) SPMUP Team 

options = struct(...
    'outdir', [BIDS_dir filesep 'derivatives'], ...
    'keep_data', 'off',  ...
    'overwrite_data', 'on', ...
    'QC', 'on', ...
    'anat',[], ...
    'task',[], ...
    'removeNvol', '0', ...
    'despike', 'on', ...
    'ignore_fieldmaps', 'off', ...
    'VDM', [], ...
    'motionexp', 'off', ...
    'scrubbing', 'on', ...
    'norm', 'T1norm', ...
    'skernel', 3, ...
    'derivatives', '1', ...
    'carpet_plot', 'on', ...
    'GLM', 'on');
