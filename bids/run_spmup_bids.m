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

%% process
try
    parpool(feature('numCores')-1); % use all available cores -1
catch no_parpool
    disp(no_parpool.message)
end

parfor s =1: numel(subjects)
    subjects{s} = spmup_BIDS_preprocess(subjects{s}, options);
    if strcmpi(options.GLM,'on')
        subjects{s} = spmup_BIDS_1rstlevel(subjects{s}, options);
    end
    disp('-----------------------')
    fprintf('subject %g finished!',s)
    disp('-----------------------')    
end
save([options.outdir filesep 'subjects.mat'],'subjects');

%% save and report QC
AnatQA = table(cellfun(@(x) x.anat_qa.SNR, subjects)',...
    cellfun(@(x) x.anat_qa.CNR, subjects)',...
    cellfun(@(x) x.anat_qa.FBER, subjects)',...
    cellfun(@(x) x.anat_qa.EFC, subjects)',...
    'VariableNames',{'SNR','CNR','FBER','EFC'});
AnatQA.Properties.Description = 'spmup AnatQC';
writetable(AnatQA,[options.outdir filesep 'AnatQC.tsv'],...
    'Delimiter','\t','FileType','text')
spmup_plotqc(AnatQA,'new')

fMRIQA = table(cellfun(@(x) x.func_qa{1}.preproc_tSNR.GM, subjects)',...
    cellfun(@(x) x.func_qa{1}.preproc_tSNR.WM, subjects)',...
    cellfun(@(x) x.func_qa{1}.preproc_tSNR.CSF, subjects)',...
    cellfun(@(x) x.func_qa{1}.preproc_tSNR.average, subjects)',...
    'VariableNames',{'tSNR GM','tSNR WM','tSNR CSF','tSNR average'});
fMRIQA.Properties.Description = 'spmup tSNR';
writetable(AnatQA,[options.outdir filesep 'fMRIQC.tsv'],...
    'Delimiter','\t','FileType','text')
spmup_plotqc(fMRIQA,'new')




