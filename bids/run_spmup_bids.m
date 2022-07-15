function subjects = run_spmup_bids(subjects, varargin)

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
% INPUT subjects is the structure coming out of spmup_BIDS_unpack
%       options is either a set of key/value pairs matching fields from 
%               spmup_getoptions or direction the options structure
%
% OUTPUT the preprocessed and analyzed data are written on disk 
%
% Example usage BIDS_dir        = 'F:\WakemanHenson_Faces\fmri';
%               options         = spmup_getoptions(BIDS_dir); 
%               [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options);
%               run_spmup_bids(subjects)
%
% Cyril Pernet July 2022
% --------------------------
%  Copyright (C) SPMUP Team 

%% Set options, matlab path and get data

% check data are there
% for s = 1:length(subjects)
%     for f = 1:size(subjects{s}.anat,1)
%         if ~exist(subjects{s}.anat{f},'file')
%             error(' %s file is missing',subjects{s}.anat{f})
%         end
%     end
%     
%     for f = 1:size(subjects{s}.func,1)
%         if ~exist(subjects{s}.func{f},'file')
%             error(' %s file is missing',subjects{s}.func{f})
%         end
%     end
%     
%     if isfield(subjects{s},'fmap')
%         for f = 1:size(subjects{s}.fmap,1)
%             if ~exist(subjects{s}.fmap{f},'file')
%             error(' %s file is missing',subjects{s}.fmap{f})
%             end
%         end
%     end
% end

% make sure all folders are in the path 
addpath(genpath(fullfile(fileparts(which('spm.m')),['toolbox' filesep 'spmup'])));

% check inputs
if nargin >= 2
    if nargin == 2 && isstruct(varargin{1})
        options = varargin{1};
    else
        for v=1:2:length(varargin)
            if any(contains(varargin{v},fieldnames(options)))
                options.(varargin{v}) = varargin{v+1};
            end
        end
    end
end

%% preprocess
try
    parpool(feature('numCores')-1); % use all available cores -1
catch no_parpool
    disp(no_parpool.message)
end

parfor s =1: numel(subjects)
    subject{s} = spmup_BIDS_preprocess(subjects{s}, options)
    if strcmpi(options.GLM,'on')
        subject{s} = spmup_BIDS_1rstlevel(subjects{s}, options)
    end
end
save subjects subjects

