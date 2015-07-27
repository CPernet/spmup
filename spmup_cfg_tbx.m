function spmup = spmup_cfg_tbx

% SPM Configuration file for the spm utility plus toolbox
% version 1 - july 2015
% ------------------------------------------------------
% Copyright (C) spmup team 2015

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','spmup')); end

% ------------------------
% generic image selection
% -----------------------
select_images         = cfg_files;
select_images.tag     = 'select_images';
select_images.name    = 'Data';
select_images.help    = {'Select image(s) to process'};
select_images.num     = [1 Inf];
select_images.filter  = 'image';
select_images.ufilter = '.*';

% ------------------------
% generic mask selection
% -----------------------
select_mask         = cfg_files;
select_mask.tag     = 'select_mask';
select_mask.name    = 'Mask image';
select_mask.help    = {'Select mask image'};
select_mask.num     = [1 1];
select_mask.filter  = 'image';
select_mask.ufilter = '.*';

% ---------------------------------------------------------------------
% auto orient menu - select image to get M from
% ---------------------------------------------------------------------
reorient_matrix         = cfg_entry;
reorient_matrix.tag     = 'reorient_matrix';
reorient_matrix.name    = 'Reorientation Matrix';
reorient_matrix.help    = {'From which image to you want to reorientation matrix from? (default = 1)'};
reorient_matrix.strtype = 'r';
reorient_matrix.num     = [1  1];

% ---------------------------------------------------------------------
% despiking menu
% ---------------------------------------------------------------------
despike_options         = cfg_branch;
despike_options.tag     = 'despike_options';
despike_options.name    = 'Options';
despike_options.val     = {select_mask};

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

% ---------------
% make branches
% ---------------
orient_br         = cfg_exbranch;
orient_br.tag     = 'autoorient';
orient_br.name    = 'Automatic Reorientation';
orient_br.val     = {select_images reorient_matrix};
orient_br.help    = {'Orient images using a simple affine registration to the template.'
    'This can be useful to set automatically the origin'};
orient_br.prog    = @spmup_reorient;
orient_br.vout    = @orient_out;



despike_br         = cfg_exbranch;
despike_br.tag     = 'despiking';
despike_br.name    = 'Despking';
despike_br.val     = {select_images despike_options};
despike_br.help    = {'Despike time series by applying a median smoother or a locally'
    'weighted smoother (lowess). Time series can be noisy with spikes'
    'and temporal smoothing can lead to better estimates'};
despike_br.prog = @spmup_despike;
despike_br.vout = @despike_out;



% ---------------------------------------------------------------------
% cfg spmup toolbox
% ---------------------------------------------------------------------
spmup          = cfg_choice;
spmup.tag      = 'spmup_cfg';
spmup.name     = 'SPM Utility Plus toolbox';
spmup.values   = {orient_br despike_br};
spmup.help     = {'Sets of utilities to get the most of your mass-univariate analyses'};

end

% -------------------------------------------------------------------------
% JOBS
% -------------------------------------------------------------------------
function M = spmup_reorient(job)
M = spmup_auto_reorient(job.select_images, job.reorient_matrix);
end

function dep = orient_out(job)
dep            = cfg_dep;
dep.sname      = sprintf('Orientation matrix from image %g',job.reorient_matrix);
end

function new_files = spmup_despike(job)
new_files = spmup_despike(job.select_images, job.select_mask, job.despike_options);
end

function dep = despike_out
dep            = cfg_dep;
dep.sname      = sprintf('Despiked images');
end

