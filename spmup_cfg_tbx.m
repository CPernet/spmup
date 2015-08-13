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

% ------------------------
% generic mat selection
% -----------------------
select_mat          = cfg_files;
select_mat.tag     = 'select_mat';
select_mat.name    = 'MAT file';
select_mat.help    = {'Select SPM.mat'};
select_mat.num     = [1 1];
select_mat.filter  = 'mat';
select_mat.ufilter = '^SPM\.mat$';

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
despike_option         = cfg_branch;
despike_option.tag     = 'despike_options';
despike_option.name    = 'Options';
despike_option.help    = {'Choose a brain mask in which voxels will be despiked'};
despike_option.val     = {select_mask};

% ---------------------------------------------------------------------
% 1st level QA
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% hrf boost
% ---------------------------------------------------------------------
boost_temporal_val         = cfg_entry;
boost_temporal_val.tag     = 'boost_temporal_val';
boost_temporal_val.name    = 'Temporal Range';
boost_temporal_val.strtype = 'r';
boost_temporal_val.num     = [1  1];
boost_temporal_val.def     = @(val)boost_default_shift; 

boost_option         = cfg_branch;
boost_option.tag     = 'boost_option';
boost_option.name    = 'Temporal range';
boost_option.help    = {'enter the default range of time dispersion allowed (default = 2 sec)'};
boost_option.val     = {boost_temporal_val};

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
despike_br.val     = {select_images despike_option};
despike_br.help    = {'Despike time series by applying a median smoother or a locally'
    'weighted smoother (lowess). Time series can be noisy with spikes'
    'and temporal smoothing can lead to better estimates'};
despike_br.prog = @spmup_despike;
despike_br.vout = @despike_out;

boost_br         = cfg_exbranch;
boost_br.tag     = 'boost';
boost_br.name    = 'Hrf boost';
boost_br.val     = {select_mat boost_option}; 
boost_br.help    = {'Hrf boost combined the beta parameters of the hrf with 1st and'
    'or the 2nd derivatives to create a boosted version which reflects the ''true'' '
    'height from the full model - for details see Pernet (2014) Front. Neurosci 8, 1 '
    '<doi: 10.3389/fnins.2014.00001> '};
boost_br.prog = @spmup_boost;


% ---------------------------------------------------------------------
% cfg spmup toolbox
% ---------------------------------------------------------------------
spmup          = cfg_choice;
spmup.tag      = 'spmup_cfg';
spmup.name     = 'SPM Utility Plus toolbox';
spmup.values   = {orient_br despike_br boost_br};
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

function dep = despike_out(job)
dep            = cfg_dep;
dep.sname      = sprintf('Despiked images');
end

function spmup_boost(job)
spmup_hrf_boost(job.select_mat, job.boost_option);
end

function val = boost_default_shift
val = 2; % seconds around the hrf model peak
end