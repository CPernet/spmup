function spmup = spmup_cfg_tbx

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
AO.name    = 'Image(s) to reorient';
AO.help    = {'Select images to realign'};
AO.num     = [1 Inf];
AO.filter  = 'image';
AO.ufilter = '.*';

% ---------------------------------------------------------------------
% despiking menu
% ---------------------------------------------------------------------
DS         = cfg_files;
DS.tag     = 'DS';
DS.name    = 'Image(s) to despike';
DS.filter = 'image';
DS.ufilter = '.*';
DS.num     = [1 Inf];
DS.help    = {'Select images to despike'};

% change the menu to get 3 things in
%       P the names of the fMRI images (time-series) - mandatory
%       M the name of the mask - optional
%       flags defines options to be used - optional; list then as defaults
%             'flags.auto_mask','off' or 'on' if M is not provided, auto_mask is
%             'on' but if set to 'off' the user is prompted to select a mask

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
orient_br.val     = {AO};
orient_br.help    = {'Orient images using a simple affine registration to the template.'
                     'This can be useful after slice timing and realignment but before coregister,'
                     'segment and normalize, applying to both structural and functional data'};                
orient_br.prog = @spmup_reorient;

despike_br         = cfg_exbranch;
despike_br.tag     = 'despiking';
despike_br.name    = 'Despking';
despike_br.val     = {DS};
despike_br.help    = {'Despike time series by applying a median smoother or a locally'
                      'weighted smoother (lowess). Time series can be noisy with spikes'
                      'and temporal smoothing can lead to better estimates'};                
despike_br.prog = @spmup_despike;

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
function spmup_reorient(job)
   spmup_auto_reorient(job);
end

function spmup_despike(job)
   spmup_despike(job);
end

