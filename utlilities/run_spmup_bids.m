% unpacks, preprocess, and run subject level GLM with the spmup toolbox

%% Set options, matlab path and get data
clear
clc

% set options
choices = struct(...
    'removeNvol', 5, ...
    'keep_data', 'off',  ...
    'overwrite_data', 'on', ...
    'despike', 'on', ...
    'drifter', 'on', ...
    'motionexp', 'on', ...
    'scrubbing', 'on', ...
    'compcor', 'on', ...
    'norm', 'T1norm', ...
    'skernel', [8 8 8], ...
    'derivatives', '1', ...
    'ignore_fieldmaps', 'on', ...
    'QC', 'off', ...
    'realign_unwarp', 'off',...
    'carpet_plot', 'off');

% add spm12 and spmup to path
addpath(genpath('D:\Dropbox\GitHub\spmup'))
addpath('D:\Dropbox\Code\MATLAB\Neuroimaging\SPM\spm12')

% to run on data sets from BIDS examples
% https://drive.google.com/drive/u/0/folders/0B2JWN60ZLkgkMGlUY3B4MXZIZW8
% BIDS_dir = 'D:\BIDS\ds001\rawdata';
BIDS_dir = 'D:\BIDS\ds003\rawdata';
% BIDS_dir = 'D:\BIDS\ds009\rawdata';

% BIDS_dir = 'D:\BIDS\7t_trt\rawdata';
% choices.acq = 'fullbrain';

% to run on local data
% BIDS_dir = 'D:\BIDS\McGurk\rawdata';
% BIDS_dir = 'D:\BIDS\AVT\rawdata';


%% copy file to derivative folder and unpack data
[BIDS,subjects,options] = spmup_BIDS_unpack(BIDS_dir,choices);


%% preprocess
for s=1:4 %numel(subjects)
% parfor s=1%:size(subjs_ls,2)
    [anatQA_grp{s}, fMRIQA, subjects, options] = ...
        spmup_BIDS_preprocess(BIDS_dir, BIDS, subjects, s, options);
    fMRIQA_grp.tSNR{s,1} = fMRIQA.tSNR;
        fMRIQA_grp.meanFD{s,1} = fMRIQA.meanFD;
end


%% run first level GLM
% for s=1:numel(subjects)
%     [] = spmup_BIDS_1rstlevel(BIDS_dir, BIDS, subjects, s, options)
% end

