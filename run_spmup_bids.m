% unpacks, preprocess, and run subject level GLM with the spmup toolbox

%% Set options, matlab path and get data
clear
clc

% set options
choices = struct(...
    'removeNvol', 1, ...
    'keep_data', 'on',  ...
    'overwrite_data', 'off', ...
    'despike', 'on', ...
    'drifter', 'off', ...
    'motionexp', 'off', ...
    'scrubbing', 'off', ...
    'compcor', 'off', ...
    'norm', 'T1norm', ...
    'skernel', [8 8 8], ...
    'derivatives', '1', ...
    'ignore_fieldmaps', 'off', ...
    'QC', 'on');
    % 'outdir', '';

% add spm12 and spmup to path
addpath(genpath('D:\Dropbox\GitHub\spmup'))
addpath('D:\Dropbox\Code\MATLAB\Neuroimaging\SPM\spm12')

% to run on data sets from BIDS examples
% https://drive.google.com/drive/u/0/folders/0B2JWN60ZLkgkMGlUY3B4MXZIZW8
BIDS_dir = 'D:\BIDS\ds001\rawdata';
% BIDS_dir = 'D:\BIDS\ds003\rawdata';
% BIDS_dir = 'D:\BIDS\ds009\rawdata';
% BIDS_dir = 'D:\BIDS\7t_trt\rawdata';

% to run on local data
% BIDS_dir = 'D:\BIDS\McGurk\rawdata';
% BIDS_dir = 'D:\BIDS\AVT\rawdata';


%% copy file to derivative folder and unpack data
[BIDS,subjects,options] = spmup_BIDS_unpack(BIDS_dir,choices);


%% preprocess
for s=1:numel(subjects)
% parfor s=1%:size(subjs_ls,2)
    spmup_BIDSjob(BIDS_dir,BIDS,subjects,s,options)
end


%% run first level GLM