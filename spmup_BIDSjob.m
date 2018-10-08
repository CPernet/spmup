function spmup_BIDSjob(BIDS_dir,BIDS,subjects,s,options,start_at)

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
% Despiking and slice timing for each run
% ----------------------------------------
for frun = 1:size(subjects{s}.func,1) % each run

    filesin = subjects{s}.func{frun};
    
    metadata = spm_BIDS(BIDS,'metadata', 'sub', subjs_ls{s}, 'run', runs_ls{frun});
    
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
    if start_at ~= 1
        
        % remove from the 4D the files we don't want and proceed
        fprintf('\nadjusting 4D file sizes run %g \n',frun)

        three_dim_files = spm_file_split(subjects{s}.func{frun});
        V = three_dim_files; V(1:start_at) = [];
        spm_file_merge(V,filesin);
        spm_unlink(three_dim_files.fname)
        
    end
    [filepath,filename,ext] = fileparts(filesin);
    
    % -----------------------------------------
    % despiking using adaptive median filter
    % -----------------------------------------
    if strcmp(options.despike,'on')
        
        if ~isfield(options,'despiking_window')
            options.despiking_window = [];
        end
        flags = struct('auto_mask','on', 'method','median', 'window', ...
            options.despiking_window,'skip',0);
        [~,filename,ext] = fileparts(filesin);
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                && ~exist([filepath filesep 'despiked_' filename ext],'file'))
            Vin = spmup_despike(filesin,[],flags);
            filesin = Vin.fname; clear Vin;
        end
    end
    
    % -------------
    % slice timing
    % -------------
    % sliceorder - vector containig the acquisition time for each slice in milliseconds
    % refslice   - time in milliseconds for the reference slice
    % timing     - [0 TR]
    
    [~,~,info_position] = intersect([filename(1:end-4) 'events.tsv'], all_names);
    
    sliceorder  = SliceTiming; % time
    refslice    = sliceorder(round(length(SliceTiming)/2));
    timing      = [0 RepetitionTime];
    if strcmp(options.despike,'on')
        st_files{frun} = [filepath filesep 'st_despiked_' filename ext]; %#ok<*AGROW>
    else
        st_files{frun} = [filepath filesep 'st_' filename ext];
    end
    
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~exist(st_files{frun},'file'))
        fprintf('\n\nstarting slice timing correction run %g subject %g \n',frun,s)
        spm_slice_timing(filesin, sliceorder, refslice, timing,'st_');
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
        [filepath,filename,ext] = fileparts(st_files{frun});
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
            matlabbatch{1}.spm.spatial.realignunwarp.data(frun).pmscan = {vdm{frun}};
        end
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality    = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep        = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm       = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm        = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp    = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap      = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight     = '';
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn   = [12 12];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda   = 100000;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm       = 0;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot      = [4 5];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot      = [];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm   = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem      = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi      = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich  = [2 1];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp  = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap     = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask     = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix   = 'ur';
        spm_jobman('run',matlabbatch); clear matlabbatch;
    end
    disp('setting up AC coordinate based on template')
    tmp = realigned_file; tmp{length(tmp)+1} = mean_realigned_file;
    spmup_auto_reorient(tmp',length(tmp)); clear tmp;
end

% ---------------
% Normalization
% ---------------

% Coregistration / Segmentation / Normalization
% ----------------------------------------------
[filepath,filename,ext] = fileparts(subjects{s}.anat);
normalization_file = [filepath filesep 'y_' filename ext];
NormalizedAnat_file  = [filepath filesep 'wm' filename ext];
class{1} = [filepath filesep 'c1' filename ext];
class{2} = [filepath filesep 'c2' filename ext];
class{3} = [filepath filesep 'c3' filename ext];
EPI_class{1} = [filepath filesep 'c1r' filename ext];
EPI_class{2} = [filepath filesep 'c2r' filename ext];
EPI_class{3} = [filepath filesep 'c3r' filename ext];
Normalized_class{1} = [filepath filesep 'wc1' filename ext];
Normalized_class{2} = [filepath filesep 'wc2' filename ext];
Normalized_class{3} = [filepath filesep 'wc3' filename ext];
for frun = 1:size(subjects{s}.func,1)
    [filepath,filename,ext] = fileparts(realigned_file{frun});
    Normalized_files{frun} = [filepath filesep 'w' filename ext];
end

if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
        && ~exist(normalization_file,'file'))
    fprintf('\n\nsubject %g: coregister, segment and normalize \n',s)
    if exist('matlabbatch','var')
        clear matlabbatch
    end
    
    % coregister anatomical to mean EPI
    % ---------------------------------
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {mean_realigned_file};
    matlabbatch{end}.spm.spatial.coreg.estimate.source = {subjects{s}.anat};
    matlabbatch{end}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
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
    
    % if field maps, normalize EPI from T1
    % -------------------------------------
    if strcmpi(options.norm,'T1norm') % ~isempty(BIDS.subjects(s).fmap)
        
        % segment the coregistered T1 (not resliced)
        % -------------------------------------------
        matlabbatch{4}.spm.spatial.preproc.channel.vols(1) = ...
        cfg_dep('Coregister: Estimate: Coregistered Images', ...
        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','cfiles'));
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
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(1) = ...
            cfg_dep('Segment: Bias Corrected (1)', ...
            substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(2) = ...
            cfg_dep('Segment: c1 Images', ...
            substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(3) = ...
            cfg_dep('Segment: c2 Images', ...
            substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(4) = ...
            cfg_dep('Segment: c3 Images', ...
            substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));

        % normalize EPI
        matlabbatch{end}.spm.spatial.normalise.write.subj(2).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        for frun = 1:size(subjects{s}.func,1)
            matlabbatch{end}.spm.spatial.normalise.write.subj(2).resample(frun,:) = {realigned_file{frun}};
        end
        
        matlabbatch{end}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70 ; 78 76 85];
        matlabbatch{end}.spm.spatial.normalise.write.woptions.vox = [2 2 2]; %might wanna change that to EPI original res
        matlabbatch{end}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{end}.spm.spatial.normalise.write.woptions.prefix = 'w';
        
    else  % normalize EPI on EPI template (old routine)
        % -------------------------------------------
        matlabbatch{4}.spm.tools.oldnorm.estwrite.subj.source(1) = {mean_realigned_file};
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
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2]; %might wanna change that to EPI original res
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.interp = 4;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
        
    end
    
    save('batch.mat', 'matlabbatch')
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
        % Basic QA for anatomical data is to get SNR, CNR, FBER and Entropy
        % This is useful to check coregistration and normalization worked fine
        tmp = spmup_anatQA(NormalizedAnat_file,Normalized_class{1},Normalized_class{2});
        save([fileparts(NormalizedAnat_file) filesep 'anatQA.mat'],'tmp');
        anatQA{s} = tmp; clear tmp
    end
    
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~exist([fileparts(Normalized_files{1}) filesep 'fMRIQA.mat'],'file'))
        fprintf('subject %g: fMRI Quality control \n',s)
        % For functional data, QA is consists in getting temporal SNR and then
        % check for motion - here we also compute additional regressors to
        % account for motion
        
        if strcmpi(options.scrubbing,'on')
            flags = struct('motion_parameters','on','globals','on','volume_distance','off','movie','off', ...
                'AC', [], 'average','on', 'T1', 'on');
        else
            flags = struct('motion_parameters','off','globals','off','volume_distance','on','movie','off', ...
                'AC', [], 'average','on', 'T1', 'on');
        end
        
        for frun = 1:size(subjects{s}.func,1)
            fMRIQA.tSNR{s, frun} = spmup_temporalSNR(stats_ready{frun},EPI_class,0);
            tmp = spmup_first_level_qa(NormalizedAnat_file,cell2mat(Normalized_files(frun)),flags);
            fMRIQA.meanFD{s,frun} = mean(spmup_FD(cell2mat(tmp))); clear tmp
            QA.tSNR = fMRIQA.tSNR{s,frun}; QA.meanFD = fMRIQA.meanFD{s,frun};
            save([fileparts(Normalized_files{frun}) filesep 'fMRIQA.mat'],'QA'); clear QA
        end
    end
end


% --------------------
% 1st level modelling
% --------------------
filepath = fileparts(fileparts([BIDS_dir filesep subjects{s}.func(frun,:)]));
SPMmat_file = [filepath filesep 'Stats' filesep 'SPM.mat'];
if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(SPMmat_file,'file'))
    
    clear matlabbatch; mkdir([filepath filesep 'Stats'])
    matlabbatch{1}.spm.stats.fmri_spec.dir = {[filepath filesep 'Stats']};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = BIDS.subjects(1).func(1).meta.RepetitionTime;
    N = length(BIDS.subjects(s).func(1).meta.SliceTiming);
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = N;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(N/ 2);
    
    % sessions, onsets and durations
    for frun = 1:size(subjects{s}.func,1)
        [~,filename,~]=fileparts(subjects{s}.func(frun,:));
        [~,~,info_position] = intersect([filename(1:end-4) 'events.tsv'], all_names);
        N_events = length(BIDS.subjects(s).func(info_position).meta.trial_type);
        for n=1:N_events
            cond{n} = cell2mat(BIDS.subjects(s).func(info_position).meta.trial_type(n,:));
        end
        cond = unique(cond);
        all_cond{frun} = cond;
        N_cond = length(cond);
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans = {stats_ready{frun}};
        onsets = NaN(N_events,N_cond);
        durations = NaN(N_events,N_cond);
        for C = 1:N_cond
            for n=1:N_events
                if strcmp(cell2mat(BIDS.subjects(s).func(info_position).meta.trial_type(n,:)),cond{C})
                    onsets(n,C) = BIDS.subjects(s).func(info_position).meta.onset(n);
                    durations(n,C) = BIDS.subjects(s).func(info_position).meta.duration(n);
                end
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).name = cond{C};
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).onset = onsets(~isnan(onsets(:,C)),C);
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).duration = durations(~isnan(durations(:,C)),C);
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).orth = 1;
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(frun).multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(frun).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(frun).multi_reg = {multi_reg{frun}};
        matlabbatch{1}.spm.stats.fmri_spec.sess(frun).hpf = 128;
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    if strcmp(options.derivatives,'off') % otherwise do it after stats
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    end
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    if BIDS.subjects(1).func(1).meta.RepetitionTime >=1
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    else
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'FAST';
    end
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', ...
        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    save([filepath filesep 'stats_batch_sub' num2str(s) '.mat'],'matlabbatch');
    spm_jobman('run',matlabbatch); clear matlabbatch;
end

% contrasts
% we can assume that at 1 contrast per condition is computed, if there
% is only one run, no need this is the same as betas, if there are
% seveal runs, simply match condition labels - and do it for each
% derivatives if any
if size(subjects{s}.func,1) > 1
    if ~exist('CI','var') %% if we did not recompute the model this is missing
        for frun = 1:size(subjects{s}.func,1)
            [~,filename,~]=fileparts(subjects{s}.func(frun,:));
            [~,~,info_position] = intersect([filename(1:end-4) 'events.tsv'], all_names);
            N_events = length(BIDS.subjects(s).func(info_position).meta.trial_type);
            for n=1:N_events
                cond{n} =cell2mat(BIDS.subjects(s).func(info_position).meta.trial_type(n,:));
            end
            all_cond{frun} = unique(cond);
            if frun > 1
                CI = intersect(all_cond{frun-1},all_cond{frun});
            end
        end
    end
    
    if ~isempty(CI)
        filepath = fileparts(SPMmat_file);
        con_file1 = [filepath filesep 'con_0001.nii'];
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                && ~exist(con_file1,'file'))
            
            % check contrast weight for run 1 and replicate if needed
            for n=1:length(BIDS.subjects(s).func(info_position).meta.trial_type)
                cond{n} = cell2mat(BIDS.subjects(s).func(info_position).meta.trial_type(n));
            end
            cond = unique(cond);
            v = zeros(length(cond));
            for c=1:length(cond) % that's the columns to span
                for a=1:length(all_cond{1}) % that's the contrast to do
                    if strcmp(all_cond{1}{a},cond{c})
                        v(a,c) = 1;
                    end
                end
            end
            
            matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_file};
            for frun = 1:size(subjects{s}.func,1)
                for c=1:length(CI)
                    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = CI{c};
                    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = v(c,:);
                    if size(BIDS.subjects(s).func,2) == 1
                        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
                    else
                        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'repl';
                    end
                end
            end
            matlabbatch{1}.spm.stats.con.delete = 1;
            spm_jobman('run',matlabbatch); clear matlabbatch;
        end
    end
    
    if strcmp(options.derivatives,'on')
       outfiles = spmup_hrf_boost(SPMmat_file);
       spmup_smooth_boostedfiles(outfiles{1},options.skernel);
       if length(outfiles) == 2
           spmup_smooth_boostedfiles(outfiles{2},options.skernel);
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

