function [anatQA, fMRIQA, subjects, options] = spmup_BIDS_preprocess(BIDS_dir, BIDS, subjects, s, options)

% routine to preprocess BIDS fMRI data - with various options available
% FORMAT spmup_BIDS_preprocess(BIDS_dir, BIDS, subjects, s)
%        spmup_BIDS_preprocess(BIDS_dir, BIDS, subjects, s, options)
%
% INPUTS 
%           - BIDS_dir is the BIDS directory
%           - BIDS: the structure returned by spm_BIDS and possibly modified
%           by spmup_BIDS_unpack
%           - subjects: a structure containing the fullpath of the unpacked anat,
%           fmap and func files for each subject (see spmup_BIDS_unpack)
%           - s is the subject index to preprocess
%           - options is a structure with the following fields:
%               .outdir = where to write the data
%               .removeNvol = number of initial volumes to remove
%               .keep_data = 'off' (default) or 'on' to keep all steps - off means
%                            only the last processed data are available
%               .overwrite_data = 'on' turning it 'off' is useful to restart
%                                  some processing while kepping previous steps
%               .QC = 'on' (default) or 'off' performs quality controls for
%                     anatomical (T1) scans and EPI time series
%               .despike = 'on' (default) or 'off' runs median despiking
%               .drifter = 'off' ('default') or 'on' removes cardiac and respiratory signals using the drifter toolbox
%               .motionexp = 'off' (default) or 'on' compute 24 motion parameters
%               .scrubbing = 'off' (default) or 'on' find outliers in motion derivatives and in globals
%               .compcor = 'on' (default) or 'off' does the equivalent of compcor
%               .norm = 'EPInorm' (default) or 'T1norm' choice of the type of template for normalization
%               .ignore_fieldmaps = 'on' or 'off' (default) to include distorsion correction for T1norm
%               .skernel = [8 8 8] by default is the smoothing kernel
%               .derivatives = 'off', 1 or 2 to use for GLM
%                              if dervatives are used, beta hrf get boosted
%                              and smoothing is performed after the GLM
%
% usage:
% choice = struct('removeNvol', 0, 'keep_data', 'off',  'overwrite_data', 'on', ...
%     'despike', 'off', 'drifter', 'off', 'motionexp', 'off', 'scrubbing', 'off', ...
%     'compcor', 'off', 'norm', 'EPInorm', 'skernel', [8 8 8], 'derivatives', 'off', ...
%     'ignore_fieldmaps', 'on',  'outdir', ['..' filesep 'derivatives' filesep 'spmup_BIDS_processed'], 'QC', 'off'); % standard SPM pipeline
% [BIDS,subjects,options]=spmup_BIDS_unpack(pwd,choice)
% for s=1:numel(subjects)
% [anatQA, fMRIQA, subjects, options] = spmup_BIDS_preprocess(BIDS_dir, BIDS, subjects, s, options)
% end

% TO DO:
% - track which task for each bold run ?
% - implement fieldmap and epi types for fieldmap modality ?
% - add an spm_check_coregistration to vizualize how the spmup_autoreorient
%   worked on the anat data?
% - change prefix name appending of preprocessed data (i.e a la SPM) to
% suffix name appending as per the BIDS derivative specs
% - create json files for each preprocessing step (see bids derivative
% specs)

if isfield(options,'removeNvol')
    start_at = options.removeNvol;
else
    start_at = 1;
end

% each subject build a job structure around matlabbatch
% this is where each subject is analysed using the various options

fprintf('\n\nrunning subject %g \n',s)
spm_root = fileparts(which('spm'));
spm_jobman('initcfg');

subjs_ls = spm_BIDS(BIDS,'subjects');
runs_ls = spm_BIDS(BIDS,'runs', 'sub', subjs_ls{s});
all_names = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
            'type', 'bold');
    
        
% ---------------------------------------        
% reorient anat file to template
% ---------------------------------------

target_dir = spm_fileparts(subjects{s}.anat);
file_exists = spm_select('FPList',target_dir,'^reorient_mat_anat.*' );
if strcmp(options.overwrite_data,'on') || ( strcmp(options.overwrite_data,'off') ...
        && isempty(file_exists) )
    
    RM = spmup_auto_reorient(subjects{s}.anat); disp(' anat reoriented'); %#ok<NASGU>
    
    % saves reorient matrix
    date_format = 'yyyy_mm_dd_HH_MM';
    saved_RM_file = fullfile(target_dir, ...
        strcat('reorient_mat_anat', datestr(now, date_format), '.mat'));
    save(saved_RM_file, 'RM')
    
    % apply reorientation to functional imges
    matlabbatch{1}.spm.util.reorient.srcfiles = {};
    for irun = 1:size(subjects{s}.func,1)
        n_vol = numel(spm_vol(subjects{s}.func{irun}));
        for i_vol = 1:n_vol
            matlabbatch{1}.spm.util.reorient.srcfiles{end+1,1} = ...
                [subjects{s}.func{irun} ',' num2str(i_vol)];
        end
    end
    matlabbatch{1}.spm.util.reorient.transform.transM = RM;
    matlabbatch{1}.spm.util.reorient.prefix = '';
    
    spm_jobman('run',matlabbatch); clear matlabbatch;
end

    
% ---------------------------------------
% Despiking and slice timing for each run
% ----------------------------------------

for frun = 1:size(subjects{s}.func,1) % each run

    filesin = subjects{s}.func{frun};
    [filepath,filename,ext] = fileparts(filesin);
    
    try
        metadata = spm_BIDS(BIDS,'metadata', 'sub', subjs_ls{s}, 'run', runs_ls{frun});
    catch
        metadata{1} = spm_BIDS(BIDS,'metadata', 'sub', subjs_ls{s}, 'type', 'bold');
    end
    hdr = spm_vol(subjects{s}.func{frun});
    epi_res = diag(hdr(1).mat);
    epi_res(end) = [];
    
    for i=1:numel(metadata)
        if isfield(metadata{i}, 'SliceTiming')
            SliceTiming = metadata{i}.SliceTiming;
        end
        if isfield(metadata{i}, 'RepetitionTime')
            RepetitionTime = metadata{i}.RepetitionTime;
        end
    end

    
    % -----------------------------------------
    % remove first volumes
    % -----------------------------------------
    
    file_exists = exist(fullfile(filepath,[filename '_removevol.json']),'file');
    if start_at ~= 1 && ( strcmp(options.overwrite_data,'on')...
            || (strcmp(options.overwrite_data,'off') && ~file_exists) )
        
        % remove from the 4D the files we don't want and proceed
        fprintf('\nadjusting 4D file sizes run %g \n',frun)

        three_dim_files = spm_file_split(filesin);
        V = three_dim_files; V(1:start_at) = [];
        spm_file_merge(V,filesin);
        spm_unlink(three_dim_files.fname)
        
        % write a json file containing the details of what volumes were
        % removed (see BIDS derivatives specs)
        spm_jsonwrite(fullfile(filepath,[filename '_removevol.json']),'')

    end
    
    
    % -----------------------------------------
    % despiking using adaptive median filter
    % -----------------------------------------
    
    if strcmp(options.despike,'on')
        
        if ~isfield(options,'despiking_window')
            options.despiking_window = [];
        end
        flags = struct('auto_mask','on', 'method','median', 'window', ...
            options.despiking_window,'skip',0);
        
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                && ~exist([filepath filesep 'despiked_' filename ext],'file'))
            Vin = spmup_despike(fullfile(filepath,[filename,ext]),[],flags);
            filesin = Vin.fname; clear Vin;
        end
        
        [filepath,filename,ext] = fileparts(filesin);
        
    end
    
    
    % -------------
    % slice timing
    % -------------
    % sliceorder - vector containig the acquisition time for each slice in milliseconds
    % refslice   - time in milliseconds for the reference slice
    % timing     - [0 TR]
    
    [~,~,info_position] = intersect([filename(1:end-4) 'events.tsv'], all_names);
    
    if exist('SliceTiming', 'var')
        sliceorder  = SliceTiming; % time
        refslice    = sliceorder(round(length(SliceTiming)/2));
        timing      = [0 RepetitionTime];
        
        st_files{frun} = [filepath filesep 'st_' filename ext]; %#ok<*AGROW>
        
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                && ~exist(st_files{frun},'file'))
            fprintf('\n\nstarting slice timing correction run %g subject %g \n',frun,s)
            spm_slice_timing(filesin, sliceorder, refslice, timing,'st_');
        end
        
    else
        st_files{frun} = fullfile(filepath, [filename ext]);
    end
       
    
    % ----------------------
    % Field Map - compute VDM
    % ----------------------
    
    if ~isempty(BIDS.subjects(s).fmap)
        fprintf('computing voxel displacement map %g subject %g \n',frun,s)
        if strcmp(options.despike,'on')
            avg = [filepath filesep 'spmup_mean_st_despiked_' filename ext];
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                    && ~exist(avg,'file'))
                spmup_bascis([filepath filesep 'despiked_st_' filename ext],'mean'); % use the mean despiked slice timed EPI for QC
            end
        else
            avg = [filepath filesep 'spmup_mean_st_' filename ext];
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                    && ~exist(avg,'file'))
                spmup_bascis([filepath filesep 'st_' filename ext],'mean'); % use the mean slice timed EPI for QC
            end
        end
        
        % output images
        if isfield(BIDS.subjects(s).fmap{frun},'phasediff')
            % two magnitudes (use only 1) and 1 phase difference image
            [~,name]=fileparts(subjects{s}.fieldmap(frun).phasediff);
            vdm{frun} = [filepath filesep 'fieldmaps' filesep 'vdm5_sc' name ext];
        else
            % two magnitide images and 2 phase images
            [~,name]=fileparts(subjects{s}.fieldmap(frun).phase1);
            vdm{frun} = [filepath filesep 'fieldmaps' filesep 'vdm5_sc' name ext];
        end
        
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                && ~exist(vdm{frun},'file'))
            % input images
            if isfield(BIDS.subjects(s).fmap{frun},'phasediff')
                % two magnitudes (use only 1) and 1 phase difference image
                matlabbatch = get_FM_workflow('phasediff');
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = ...
                    {[BIDS_dir filesep subjects{s}.fieldmap(frun).phasediff]};
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = ...
                    {[BIDS_dir filesep subjects{s}.fieldmap(frun).mag1]};
            else
                % two magnitide images and 2 phase images
                matlabbatch = get_FM_workflow('phase&mag');
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase = ...
                    {[BIDS_dir filesep subjects{s}.fieldmap.phase1]};
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag   = ...
                    {[BIDS_dir filesep subjects{s}.fieldmap.mag1]};
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase  = ...
                    {[BIDS_dir filesep subjects{s}.fieldmap.phase2]};
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag    = ...
                    {[BIDS_dir filesep subjects{s}.fieldmap.mag2]};
            end
            
            % update parameters
            echotimes =  1000.*[BIDS.subjects(s).fmap{frun}.meta.EchoTime1 BIDS.subjects(s).fmap{frun}.meta.EchoTime2]; % change to ms
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = echotimes;
            
            if isfield(BIDS.subjects(s).fmap{frun}.meta,'TotalReadoutTime')
                TotalReadoutTime = BIDS.subjects(s).fmap{frun}.meta.TotalReadoutTime;
            elseif isfield(BIDS.subjects(s).fmap{frun}.meta,'RepetitionTime')
                TotalReadoutTime = BIDS.subjects(s).fmap{frun}.meta.RepetitionTime;
            elseif isfield(BIDS.subjects(s).fmap{frun}.meta,'EffectiveEchoSpacing')
                TotalReadoutTime = (BIDS.subjects(s).func{frun}.meta.NumberOfEchos-1)*BIDS.subjects(s).func{frun}.meta.EffectiveEchoSpacing;
            end
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = TotalReadoutTime;
            
            if isfield(BIDS.subjects(s).fmap{frun}.meta,'PulseSequenceType')
                if sum(findstr(BIDS.subjects(s).fmap{frun}.meta.PulseSequenceType,'EPI')) ~= 0
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
                else
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
                end
            else
                disp('using default sequence! assuming non-EPI acquisition')
            end
            
            if isfield(BIDS.subjects(s).fmap{frun}.meta,'PhaseEncodingDirection')
                if strcmp(BIDS.subjects(s).fmap{frun}.meta.PhaseEncodingDirection,'j') ...
                        || strcmp(BIDS.subjects(s).fmap{frun}.meta.PhaseEncodingDirection,'y')
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = 1;
                elseif strcmp(BIDS.subjects(s).fmap{frun}.meta.PhaseEncodingDirection,'-j') ...
                        || strcmp(BIDS.subjects(s).fmap{frun}.meta.PhaseEncodingDirection,'-y')
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = -1;
                end
            elseif isfield(BIDS.subjects(s).func(frun).meta,'PhaseEncodingDirection')
                if strcmp(BIDS.subjects(s).func(frun).meta.PhaseEncodingDirection,'j') ...
                        || strcmp(BIDS.subjects(s).func(frun).meta.PhaseEncodingDirection,'y')
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = 1;
                elseif strcmp(BIDS.subjects(s).func(frun).meta.PhaseEncodingDirection,'-j') ...
                        || strcmp(BIDS.subjects(s).func(frun).meta.PhaseEncodingDirection,'-y')
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = -1;
                end
            else
                error('No phase encoding direction found')
            end
            
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {avg}; % use the mean despiked / slice timed image
            spm_jobman('run',matlabbatch); clear matlabbatch;
        end
    end
    
    % cleanup
    if strcmpi(options.keep_data,'off') 
        delete(filesin); % original or despiked
        if strcmp(options.despike,'on')
            [filepath,filename,ext]=fileparts(filesin);
            delete([filepath filesep filename(10:end) ext]);
        end
    end
    
end % end processing per run


% -----------------------
% Realignment across runs
% ------------------------

if isempty(BIDS.subjects(s).fmap)
    
    % file out of realign will be
    [filepath,filename,ext] = fileparts(st_files{1});
    mean_realigned_file = [filepath filesep 'mean' filename ext];
    for frun = 1:size(subjects{s}.func,1)
        realigned_file{frun} = st_files{frun}; % because we don't reslice, simple encode the linear transform in the header
        [filepath,filename] = fileparts(st_files{frun});
        multi_reg{frun} = [filepath filesep 'rp_' filename '.txt'];
    end
    
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~exist(mean_realigned_file,'file'))
        fprintf('subject %g: starting realignment \n',s)
        for frun = 1:size(subjects{s}.func,1)
            matlabbatch{1}.spm.spatial.realign.estwrite.data{frun} = {st_files{frun}};
        end
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.quality = 1;
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.interp  = 4;
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.weight  = '';
        matlabbatch{end}.spm.spatial.realign.estwrite.roptions.which   = [0 1]; %only reslice the mean
        matlabbatch{end}.spm.spatial.realign.estwrite.roptions.interp  = 4;
        matlabbatch{end}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
        matlabbatch{end}.spm.spatial.realign.estwrite.roptions.mask    = 1;
        matlabbatch{end}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';
        
        spm_jobman('run',matlabbatch); clear matlabbatch;
        
    end
else
    fprintf('subject %g: starting realignment and unwarping \n',s)
    % file out of realign will be
    [~,filename,ext] = fileparts(st_files{1});
    mean_realigned_file = [filepath filesep 'meanur' filename ext];
    for frun = 1:size(subjects{s}.func,1)
        [~,filename,ext] = fileparts(st_files{frun});
        realigned_file{frun} = [filepath filesep 'ur' filename ext]; % because we have the reslice here (not linear)
        multi_reg{frun} = [filepath filesep 'rp_' filename '.txt'];
    end
    
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~exist(mean_realigned_file,'file'))
        for frun = 1:size(subjects{s}.func,1)
            matlabbatch{1}.spm.spatial.realignunwarp.data(frun).scans = {st_files{frun}};
            matlabbatch{end}.spm.spatial.realignunwarp.data(frun).pmscan = {vdm{frun}};
        end
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.quality    = 1;
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.sep        = 4;
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.fwhm       = 5;
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.rtm        = 1;
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.einterp    = 4;
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.ewrap      = [0 0 0];
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.weight     = '';
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.basfcn   = [12 12];
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.lambda   = 100000;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.jm       = 0;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.fot      = [4 5];
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.sot      = [];
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.uwfwhm   = 4;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.rem      = 1;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.noi      = 5;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{end}.spm.spatial.realignunwarp.uwroptions.uwwhich  = [2 1];
        matlabbatch{end}.spm.spatial.realignunwarp.uwroptions.rinterp  = 4;
        matlabbatch{end}.spm.spatial.realignunwarp.uwroptions.wrap     = [0 0 0];
        matlabbatch{end}.spm.spatial.realignunwarp.uwroptions.mask     = 1;
        matlabbatch{end}.spm.spatial.realignunwarp.uwroptions.prefix   = 'ur';
        
        spm_jobman('run',matlabbatch); clear matlabbatch;

    end

end


% ---------------
% Coregistration / Segmentation
% ---------------

[filepath,filename,ext] = fileparts(subjects{s}.anat);
EPI_class{1} = [filepath filesep 'c1r' filename ext];
EPI_class{2} = [filepath filesep 'c2r' filename ext];
EPI_class{3} = [filepath filesep 'c3r' filename ext];
if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
        && ~exist(EPI_class{1},'file'))
    fprintf('\n\nsubject %g: coregister, segment \n',s)
    if exist('matlabbatch','var')
        clear matlabbatch
    end

    % coregister anatomical to mean EPI
    % ---------------------------------
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {mean_realigned_file};
    matlabbatch{end}.spm.spatial.coreg.estimate.source = {subjects{s}.anat};
    matlabbatch{end}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep = [16 8 4 2 1];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol = ...
        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    % reslice anatomical to mean EPI
    matlabbatch{2}.spm.spatial.coreg.write.ref = {mean_realigned_file};
    matlabbatch{end}.spm.spatial.coreg.write.source = {subjects{s}.anat};
    matlabbatch{end}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{end}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{end}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{end}.spm.spatial.coreg.write.roptions.prefix = 'r';
    
    
    % run the segmentation on the resliced T1 to get tissue classes
    % in the same space as the EPI data before mormalization
    % -------------------------------------------------------
    matlabbatch{3}.spm.spatial.preproc.channel.vols(1) = ...
        cfg_dep('Coregister: Reslice: Resliced Images', ...
        substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','rfiles'));  
    matlabbatch{end}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{end}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{end}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(1).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,1']};
    matlabbatch{end}.spm.spatial.preproc.tissue(1).ngaus = 2;
    matlabbatch{end}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(2).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,2']};
    matlabbatch{end}.spm.spatial.preproc.tissue(2).ngaus = 2;
    matlabbatch{end}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(3).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,3']};
    matlabbatch{end}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{end}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(4).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,4']};
    matlabbatch{end}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{end}.spm.spatial.preproc.tissue(4).native = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(5).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,5']};
    matlabbatch{end}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{end}.spm.spatial.preproc.tissue(5).native = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(6).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,6']};
    matlabbatch{end}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{end}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{end}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{end}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{end}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{end}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{end}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{end}.spm.spatial.preproc.warp.write = [0 0];
    
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
end


% ---------------
% Normalization
% ---------------

[filepath,filename,ext] = fileparts(subjects{s}.anat);
% Normalization_file = [filepath filesep 'y_' filename ext];
NormalizedAnat_file  = [filepath filesep 'wm' filename ext];
% NormalizedReslicedAnat_file  = [filepath filesep 'wr' filename ext];
class{1} = [filepath filesep 'c1' filename ext];
class{2} = [filepath filesep 'c2' filename ext];
class{3} = [filepath filesep 'c3' filename ext];
Normalized_class{1} = [filepath filesep 'wc1' filename ext];
Normalized_class{2} = [filepath filesep 'wc2' filename ext];
Normalized_class{3} = [filepath filesep 'wc3' filename ext];
for frun = 1:size(subjects{s}.func,1)
    [filepath,filename,ext] = fileparts(realigned_file{frun});
    Normalized_files{frun} = [filepath filesep 'w' filename ext];
end

if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
        && ~exist(Normalized_files{end},'file'))
    fprintf('\n\nsubject %g: normalize \n',s)
    if exist('matlabbatch','var')
        clear matlabbatch
    end
    
    % if field maps, normalize EPI from T1
    % -------------------------------------
    if strcmpi(options.norm,'T1norm') % ~isempty(BIDS.subjects(s).fmap)
        
        % segment the coregistered T1 (not resliced)
        % -------------------------------------------
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {subjects{s}.anat};
        matlabbatch{end}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{end}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{end}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{end}.spm.spatial.preproc.tissue(1).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,1']};
        matlabbatch{end}.spm.spatial.preproc.tissue(1).ngaus = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(2).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,2']};
        matlabbatch{end}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(3).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,3']};
        matlabbatch{end}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(4).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,4']};
        matlabbatch{end}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{end}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(5).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,5']};
        matlabbatch{end}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{end}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(6).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,6']};
        matlabbatch{end}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{end}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{end}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{end}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{end}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{end}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{end}.spm.spatial.preproc.warp.write = [1 1];
        
        % normalize EPI using T1 info
        % ----------------------------
        % normalize T1 
        matlabbatch{2}.spm.spatial.normalise.write.subj(1).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(1) = ...
            cfg_dep('Segment: Bias Corrected (1)', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(2) = ...
            cfg_dep('Segment: c1 Images', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(3) = ...
            cfg_dep('Segment: c2 Images', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(4) = ...
            cfg_dep('Segment: c3 Images', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));

        % normalize EPI
        matlabbatch{end}.spm.spatial.normalise.write.subj(2).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        for frun = 1:size(subjects{s}.func,1)
            matlabbatch{end}.spm.spatial.normalise.write.subj(2).resample(frun,:) = {realigned_file{frun}};
        end
        matlabbatch{end}.spm.spatial.normalise.write.subj(2).resample(end+1,:) = {mean_realigned_file}; % adding mean image
        
        matlabbatch{end}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70 ; 78 76 85];
        matlabbatch{end}.spm.spatial.normalise.write.woptions.vox = epi_res;
        matlabbatch{end}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{end}.spm.spatial.normalise.write.woptions.prefix = 'w';
        
    else  % normalize EPI on EPI template (old routine)
        % -------------------------------------------
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source(1) = {mean_realigned_file};
        matlabbatch{end}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
        for frun = 1:size(subjects{s}.func,1)
            matlabbatch{end}.spm.tools.oldnorm.estwrite.subj.resample(frun,:) = {realigned_file{frun}};
        end
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.template = ...
            {[spm_root filesep 'toolbox' filesep 'OldNorm' filesep 'EPI.nii,1']};
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -70 ; 78 76 85];
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.vox = epi_res;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.interp = 4;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
        
    end
%     save('batch.mat', 'matlabbatch')
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
end

if strcmpi(options.keep_data,'off')
    for frun = size(subjects{s}.func,1):-1:1
        delete(realigned_file{frun});
    end
end


% --------------
% Smoothing
% --------------

if strcmp(options.derivatives,'off') % otherwise do it after stats
    for frun = 1:size(subjects{s}.func,1)
        [filepath,filename,ext] = fileparts(Normalized_files{frun});
        stats_ready{frun} = [filepath filesep 's' filename ext];
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                && ~exist(stats_ready{frun},'file'))
            fprintf('subject %g: smoothing run %g \n',s,frun);
            spm_smooth(Normalized_files{frun},stats_ready{frun},options.skernel);
        end
    end
    
    if strcmpi(options.keep_data,'off')
        for frun = size(subjects{s}.func,1):-1:1
            delete(Normalized_files{frun});
        end
    end
    
else
    stats_ready = Normalized_files;
end


% ----------------------------
% QC and additional regressors
% ----------------------------

if strcmp(options.QC,'on') %
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~exist([fileparts(NormalizedAnat_file) filesep 'anatQA.mat'],'file'))
        fprintf('subject %g: Anatomical Quality control \n',s)
        % sanity check that all images are in the same space.
        V_to_check = Normalized_class';
        V_to_check{end+1} = NormalizedAnat_file;
        spm_check_orientations(spm_vol(char(V_to_check)));
        % Basic QA for anatomical data is to get SNR, CNR, FBER and Entropy
        % This is useful to check coregistration and normalization worked fine
        tmp = spmup_anatQA(NormalizedAnat_file,Normalized_class{1},Normalized_class{2});
        save([fileparts(NormalizedAnat_file) filesep 'anatQA.mat'],'tmp');
        anatQA{s} = tmp; clear tmp
    end
    
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~exist([fileparts(Normalized_files{end}) filesep 'fMRIQA.mat'],'file'))
        fprintf('subject %g: fMRI Quality control \n',s)
        % For functional data, QA is consists in getting temporal SNR and then
        % check for motion - here we also compute additional regressors to
        % account for motion
        
        davg = spmup_comp_dist2surf(subjects{s}.anat);
        
        if strcmpi(options.scrubbing,'on')
            flags = struct('motion_parameters','on','globals','on','volume_distance','off','movie','off', ...
                'AC', [], 'average','on', 'T1', 'on');
        else
            flags = struct('motion_parameters','off','globals','off','volume_distance','on','movie','off', ...
                'AC', [], 'average','on', 'T1', 'on');
        end
        
        for frun = 1:size(subjects{s}.func,1)
            
            % sanity check that all images are in the same space.
            V_to_check = Normalized_class';
            V_to_check{end+1} = stats_ready{frun};
            spm_check_orientations(spm_vol(char(V_to_check)));
            
            fMRIQA.tSNR{s, frun} = spmup_temporalSNR(Normalized_files{frun},vNormalized_class,v0);
            
            tmp = spmup_first_level_qa(NormalizedAnat_file,vcell2mat(stats_ready(frun)),vflags);  
            fMRIQA.meanFD{s,frun} = mean(spmup_FD(cell2mat(tmp),vdavg)); 
            clear tmp
            
            QA.tSNR = fMRIQA.tSNR{s,vfrun}; 
            QA.meanFD = fMRIQA.meanFD{s,vfrun};
            save([fileparts(Normalized_files{frun}) filesep 'fMRIQA.mat'],'QA'); 
            clear QA
            
            fprintf('subject %g: fMRI Quality control: carpet plot \n',s)
            P = subjects{s}.func{frun};
            c1 = EPI_class{1};
            c2 = EPI_class{2};
            c3 = EPI_class{3};
            spmup_timeseriesplot(P, c1, c2, c3, 'motion','on','nuisances','on','correlation','on');
            
        end
    end
end


end

%% field map sub-routine
function matlabbatch = get_FM_workflow(type)

if strcmpi(type,'phasediff')
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.presubphasemag.phase = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.presubphasemag.magnitude = '';
elseif strcmpi(type,'phase&mag')
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.phasemag.shortphase = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.phasemag.shortmag = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.phasemag.longphase = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).data.phasemag.longmag = '';
else
    error('unkmown input type for field map workflow')
end

FM_template = fullfile(spm_file(which('spm'),'fpath'),'toolbox','FieldMap','T1.nii');
UF = struct('method','Mark3D','fwhm',10,'pad',0,'ws',1);
MF = struct('template','','fwhm',5,'nerode',2,'ndilate',4,'thresh',0.5,'reg',0.02);
defaultsval = struct('et',[NaN NaN],'maskbrain',1,'blipdir',1,'tert','','epifm',0,'ajm',0,'uflags',UF,'mflags',MF);
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval = defaultsval;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {FM_template};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).session.epi = '';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).matchvdm = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).sessname = 'run';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).writeunwarped = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).anat = '';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).matchanat = 0;

end

