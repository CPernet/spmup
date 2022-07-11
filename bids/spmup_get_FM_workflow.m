function matlabbatch = spmup_get_FM_workflow(type)

% routine to get the SPM batch workflow for FieldMaps
% 
% FORMAT matlabbatch = get_FM_workflow(type)
%
% INPUT type can be 'phasediff' or 'phase&mag'
% OUTPUT matlabbatch is a structure pre-filled with SPM options
%
% Cyril Pernet & Remi Gau
% --------------------------
%  Copyright (C) SPMUP Team 


if strcmpi(type,'phasediff')
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.presubphasemag.phase     = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.presubphasemag.magnitude = '';
elseif strcmpi(type,'phase&mag')
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.phasemag.shortphase      = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.phasemag.shortmag        = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.phasemag.longphase       = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.phasemag.longmag         = '';
else
    error('unkmown input type for field map workflow')
end

FM_template = fullfile(spm_file(which('spm'),'fpath'),'toolbox','FieldMap','T1.nii');
UF          = struct('method','Mark3D','fwhm',10,'pad',0,'ws',1);
MF          = struct('template','','fwhm',5,'nerode',2,'ndilate',4,'thresh',0.5,'reg',0.02);
defaultsval = struct('et',[NaN NaN],'maskbrain',1,'blipdir',1,'tert','','epifm',0,'ajm',0,'uflags',UF,'mflags',MF);

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval              = defaultsval;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {FM_template};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).session.epi                       = '';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).matchvdm                          = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).sessname                          = 'run';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).writeunwarped                     = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).anat                              = '';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).matchanat                         = 0;

