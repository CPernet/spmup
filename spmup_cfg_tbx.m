function spmup_jobs = spmup_cfg_tbx

% SPM Configuration file for the spm utility plus toolbox
% version 1 - july 2015
% ------------------------------------------------------
% Copyright (C) spmup team 2015

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','spmup')); end

% ---------------------------------------------------------------------
% auto orient menu
% ---------------------------------------------------------------------
AO         = cfg_files;
AO.tag     = 'AO';
AO.name    = 'Automatic Reorientation';
AO.filter = 'image';
AO.ufilter = '.*';
AO.num     = [1 Inf];
AO.help    = {'Align images with the template and reset (0,0,0) location'};

% ---------------------------------------------------------------------
% despiking menu
% ---------------------------------------------------------------------
DS         = cfg_files;
DS.tag     = 'DS';
DS.name    = 'Despking';
DS.filter = 'image';
DS.ufilter = '.*';
DS.num     = [1 Inf];
DS.help    = {'Despike time series by applying a median smoother'};

% ---------------------------------------------------------------------
% 1st level QA
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% hrf boost
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% percentage signal change
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% 2nd level QA
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% cfg spmup toolbox
% ---------------------------------------------------------------------
spmup_jobs          = cfg_exbranch;
spmup_jobs.tag      = 'spmup_cfg';
spmup_jobs.name     = 'SPM Utility Plus toolbox';
spmup_jobs.val      = {AO DS};
spmup_jobs.help     = {'Sets of utilities to get the most of your mass-univariate analyses'};
% spmup_jobs.prog     = @spmup_run;



