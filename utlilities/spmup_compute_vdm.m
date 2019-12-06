function vdm_files = spmup_compute_vdm(fmap_folder)

% create vdm maps for each fmri session - assumes data follow BIDS 
% 
% FORMAT: vdm_files = NHS_compute_vdm(fmap_folder)
%
% INPUT: fmap_folder is the folder where to find gradient field maps
%
% OUTPUT: vdm_files is a cell array with the names for the vdm maps for each fmri
% sessions
%
% Note: this works for Siemens - not tested on other platforms
%
% Cyril Pernet - Univertsity of Edinburgh, August 2019
% ----------------------------------------------------

spm_jobman('initcfg');
spm('Defaults','fmri')
spm_root = fileparts(which('spm.m'));

% get parameters from json files
% --------------------------------
param_files = dir([fmap_folder filesep '*.json']); 
for p=1:3
    if contains(param_files(p).name,'phase')
        fieldmap_param = spm_jsonread([param_files(p).folder filesep param_files(p).name]);
        if strcmp(fieldmap_param.PhaseEncodingDirection,'i')
            blip = 1;
        else
            blip = -1;
        end
        phase_img = [param_files(p).folder filesep param_files(p).name(1:end-5) '.nii'];
    elseif contains(param_files(p).name,'magnitude1')
        magnitude_img = [param_files(p).folder filesep param_files(p).name(1:end-5) '.nii'];
        echo1 = spm_jsonread([param_files(p).folder filesep param_files(p).name]);
        echo1 = echo1.EchoTime * 1000;
    elseif contains(param_files(p).name,'magnitude2')
        echo2 = spm_jsonread([param_files(p).folder filesep param_files(p).name]);
        echo2 = echo2.EchoTime * 1000;
    end
end

% make the batch
% --------------
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase(1)         = {phase_img};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude(1)     = {magnitude_img};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et              = [echo1 echo2];
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain       = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir         = blip;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert            = 1/fieldmap_param.PixelBandwidth*1000;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm           = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm             = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method   = 'Mark3D';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm     = 10;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad      = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws       = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {[spm_root filesep 'toolbox' filesep 'FieldMap' filesep 'T1.nii']};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm     = 5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode   = 2;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate  = 4;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh   = 0.5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg      = 0.02;

func_folder = [fileparts(fmap_folder) filesep 'func'];
time_series = dir([func_folder filesep '*nii']);
for ts = 1:size(time_series,1)
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(ts).epi = {[time_series(ts).folder filesep time_series(ts).name ',4']};
end
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm                             = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname                             = 'session';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped                        = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat                                 = '';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat                            = 0;

%% run the batch and get data out
% ---------------------------------
spm_jobman('run',matlabbatch);

% graphic for QA
if exist(['spm_' datestr(now,'yyyymmmdd') '.ps'],'file')
    movefile(['spm_' datestr(now,'yyyymmmdd') '.ps'],fullfile(fileparts(fmap_folder),'derivatives','func','FieldMapsQA.ps'));   
end

% move all created files to derivatives
vdm_files = {};
if ~contains(fmap_folder,'derivatives')
   mkdir([fileparts(fmap_folder) filesep 'derivatives' filesep 'fmap'])
   all_data = dir(fmap_folder); move_index = 1;
   for ad = 3:size(all_data,1)
       if ~strcmp(all_data(ad).name(1:4),'sub-')
           vdm_files{move_index} = fullfile(fileparts(fmap_folder),'derivatives','fmap',all_data(ad).name);
           movefile(fullfile(all_data(ad).folder, all_data(ad).name),vdm_files{move_index});
           move_index = move_index + 1;
       end
   end
end

% cleanup
func_data = dir(func_folder);
for fd = 3:size(func_data,1)
    if ~strcmp(func_data(fd).name(1:4),'sub-')
        delete(fullfile(func_data(fd).folder,func_data(fd).name));
    end
end

