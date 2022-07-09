function run_spmup_bids(BIDS,subjects, varargin)

% SPMUP can be run using SPM batch mode but you can also analyze an entire
% fMRI BIDS dataset automatically using this function - it will unpack zip
% files, preprocess (despike, slice time, realign, smooth), QA (anat, fMRI),
% save the augmented regressors for motion, censoring, etc, and run subject 
% level GLM (augmented regressors, amplitude correction, delay estimation )
% using almost all the gizmos of the toolbox. Note the GLM will be the
% simplest possible version, i.e. 1 regressor per unique events, which
% might not be suitable, although using contrasts any other combinations 
% can, in theory, be obtained.
%
% FORMAT run_spmup_bids(BIDS_dir, options)
%
% INPUT BIDS_dir is the BIDS directory to analyse
%       options are a set of key/value pairs updating the defaults
%               see spmup_getoptions
%
% OUTPUT the preprocessed and analzed data are written on disk 
%
% Example usage [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options);
%               run_spmup_bids(BIDS,subjects, 'GLM', 'off')
%
% Cyril Pernet July 2022
% --------------------------
%  Copyright (C) SPMUP Team 

%% Set options, matlab path and get data

% check data are there
for s = 1:length(subjects)
    for f = 1:size(subjects{s}.anat,1)
        if ~exist(subjects{s}.anat{f},'file')
            error(' %s file is missing',subjects{s}.anat{f})
        end
    end
    
    for f = 1:size(subjects{s}.func,1)
        if ~exist(subjects{s}.func{f},'file')
            error(' %s file is missing',subjects{s}.func{f})
        end
    end
    
    if isfield(subjects{s},'fmap')
        for f = 1:size(subjects{s}.fmap,1)
            if ~exist(subjects{s}.fmap{f},'file')
            error(' %s file is missing',subjects{s}.fmap{f})
            end
        end
    end
end

% make sure all folders are in the path 
addpath(genpath(fullfile(fileparts(which('spm.m')),['toolbox' filesep 'spmup'])));

% set default options
options = spmup_getoptions(BIDS_dir);

% check inputs
if nargin > 2
    for v=1:2:length(varargin)
        if any(contains(varargin{v},fieldnames(options)))
            options.(varargin{v}) = varargin{v+1};
        end
    end
end

if strcmpi(options.norm,'T1norm') && isempty(options.VDM)
    warning(' Calhoun et al. (2017) showed T1 normalization is not good without FieldMaps \n (%s) add VDM map to options \n swithing to ''EPInorm''', ...
        'https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.23737')
    options.norm = 'EPI';
end

%% preprocess
try
    parpool(feature('numCores')-1); % use all available cores -1
catch no_parpool
    disp(no_parpool.message)
end

anatQA_grp = cell(numel(subjects),1);
fMRIQA_grp = cell(numel(subjects),1);
updated_subjects = cell(numel(subjects),1);
parfor s =1:numel(subjects)
    [anatQA_grp{s}, fMRIQA_grp{s}, updated_subjects{s}] = ...
        spmup_BIDS_preprocess(BIDS, subjects, s, options);
end
subjects = updated_subjects;
clear updated_subjects


%% run first level GLM
% for s=1:numel(subjects)
%     [] = spmup_BIDS_1rstlevel(BIDS_dir, BIDS, subjects, s, options)
% end

