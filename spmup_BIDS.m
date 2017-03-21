function [anatQA,fMRIQA]=spmup_BIDS(BIDS_dir,choices)

% routine to read and unpack BIDS fMRI and preprocess data, build the first
% level analysis - with various options available
%
% FORMAT spmup_BIDS
%        spmup_BIDS(BIDS_dir)
%        spmup_BIDS(BIDS_dir,options)
%
% INPUTS BIDS_dir is the BIDS directory
%        choices a structure with the following fields:
%               .outdir = where to write the data
%               .removeNvol = number of initial volumes to remove
%               .keep_data = 'off' (default) or 'on' to keep all steps - off means
%                            only the last processed data are available
%               .overwrite_data = 'on' turning it 'off; is useful to restart
%                                  some processing while kepping previous steps
%               .QC = 'on' (default) or 'off' performs quality controls for
%                     anatomical (T1) scans and EPI time series
%               .despike = 'on' (default) or 'off' runs median despiking
%               .drifter = 'off' ('default') or 'on' removes cardiac and respiratory signals using the drifter toolbox
%               .motionexp = 'off' (default) or 'on' compute 24 motion parameters
%               .scrubbing = 'off' (default) or 'on' find outliers in motion derivatives and in globals
%               .compcor = 'on' (default) or 'off' does the equivalent of compcor
%               .skernel = [6 6 6] by default is the smoothing kernel
%               .derivatives = 'off', 1 or 2 to use for GLM
%                              if dervatives are used, beta hrf get boosted
%                              and smoothing is performed after the GLM
%
% usage:
% options = struct('removeNvol', 0, 'keep_data', 'on',  'overwrite_data', 'on', ...
%     'despike', 'off', 'drifter', 'off', 'motionexp', 'off', 'scrubbing', 'off', ...
%     'compcor', 'off', 'skernel', [6 6 6], 'derivatives', 1, 'ignore_fieldmaps', 'on', ...
%     'outdir', 'spmup_BIDS_processed_basic', 'QC', 'off'); % standard SPM pipeline
% [anatQA,fMRIQA]=spmup_BIDS(pwd,options)
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox


if nargin == 0
    BIDS_dir= uigetdir(pwd,'select BIDS directory');
    if BIDS_dir == 0
        return
    end
end

options = get_all_options(BIDS_dir);
I=intersect(fieldnames(options),fieldnames(choices));
for f=1:length(I)
    options = setfield(options,I{f},getfield(choices,I{f}));
end

% if ~isfield(options,'removeNvol') || isempty(options.removeNvol)
%     options.removeNvol = input('How many initial EPI volumes to discard?: '); disp(' ')
% end

spm('defaults', 'FMRI');

%% get the data info
% -------------------------------------------------------------------------
cd(BIDS_dir)
disp('getting BIDS info')
BIDS = spm_BIDS(BIDS_dir);
if isfield(options,'ignore_fieldmaps')
    if strcmp(options.ignore_fieldmaps, 'on')
        for s=1:length(BIDS.subjects)
            BIDS.subjects(s).fmap = [];
        end
    end
end

anatQA = [];
fMRIQA = [];

%% unpack data
% -------------------------------------------------------------------------
if ~isfield(options,'outdir')
    options.outdir = [BIDS_dir filesep 'spmup_BIDS_processed'];
end

mkdir(options.outdir);
% this needs update for longitinal dataset ie multiple 'sessions'
for s=1:size(BIDS.subjects,2)
    fprintf('subject %g: unpacking anatomical data \n',s)
    in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat.filename];
    subjects{s}.anat = [options.outdir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat.filename(1:end-3)];
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.anat,'file'))
        gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'anat' ]);
        spmup_auto_reorient(subjects{s}.anat); disp('anat reoriented');
    end
end

if size(BIDS.subjects(s).func,2) >= 2
    try
        parpool(feature('numCores')-1); % use all available cores -1
    catch no_parpool
        disp(no_parpool.message)
    end
end

parfor s=1:size(BIDS.subjects,2)
    if ~isempty(BIDS.subjects(s).fmap)
        fprintf('subject %g: unpacking functional data and field maps \n',s)
    else
        fprintf('subject %g: unpacking functional data \n',s)
    end
    
    for frun = 1:size(BIDS.subjects(s).func,2)
        % functional
        in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'func' filesep BIDS.subjects(s).func(frun).filename];
        subjects{s}.func(frun,:) = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep BIDS.subjects(s).func(frun).filename(1:end-3)];
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.func(frun,:),'file'))
            gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun)]);
        end
        
        if ~isempty(BIDS.subjects(s).fmap)
            % field maps
            in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap{frun}.magnitude1];
            subjects{s}.fieldmap(frun,:).mag1 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.magnitude1(1:end-3)];
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).mag1,'file'))
                gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
            end
            in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap{frun}.magnitude2];
            subjects{s}.fieldmap(frun,:).mag2 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.magnitude2(1:end-3)];
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).mag2,'file'))
                gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
            end
            
            if isfield(BIDS.subjects(s).fmap{frun},'phasediff')
                in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap{frun}.phasediff];
                subjects{s}.fieldmap(frun,:).phasediff = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.phasediff(1:end-3)];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).phasediff,'file'))
                    gunzip(in,[options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                end
            else
                in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap{frun}.phase1];
                subjects{s}.fieldmap(frun,:).phase1 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.phase1(1:end-3)];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).phase1,'file'))
                    gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                end
                in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap{frun}.phase2];
                subjects{s}.fieldmap(frun,:).phase2 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.phase2(1:end-3)];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).phase2,'file'))
                    gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                end
            end
        end
    end
end

try delete(gcp('nocreate')); end
disp('spmup has finished unpacking data')


%% run preprocessing using options
% -------------------------------------------------------------------------
spm_jobman('initcfg');
spm_root = fileparts(which('spm'));

% do the computation if needed
for s=1:size(BIDS.subjects,2)
    % each subject build a job structure around matlabbatch
    
    for frun = 1:size(BIDS.subjects(s).func,2) % each run
        filesin = [BIDS_dir filesep subjects{s}.func(frun,:)];
        [filepath,filename,ext] = fileparts(filesin);
        
        % -----------------------------------------
        % despiking using adaptive median filter
        % -----------------------------------------
        if strcmp(options.despike,'on')
            
            if ~isfield(options,'despiking_window')
                options.despiking_window = [];
            end
            flags = struct('auto_mask','on', 'method','median', 'window', ...
                options.despiking_window,'skip',options.removeNvol);
            filesin = [filepath filesep 'despiked_' filename ext];
            [~,filename,ext] = fileparts(filesin);
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(filesin,'file'))
                spmup_despike(filesin,[],flags);
            end
            
            if strcmpi(options.keep_data,'off')
                delete(filesin);
            end
        end
        
        % -------------
        % slice timing
        % -------------
        % sliceorder - vector containig the acquisition time for each slice in milliseconds
        % refslice   - time in milliseconds for the reference slice
        % timing     - [0 TR]
        
        sliceorder  = BIDS.subjects(s).func(frun).meta.SliceTiming; % time
        refslice    = sliceorder(round(length(BIDS.subjects(s).func(frun).meta.SliceTiming)/2)); % time
        timing      = [0 BIDS.subjects(s).func(frun).meta.RepetitionTime];
        if strcmp(options.despike,'on')
            st_files{frun} = [filepath filesep 'st_despiked_' filename ext];
        else
            st_files{frun} = [filepath filesep 'st_' filename ext];
        end
        
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(st_files{frun},'file'))
            fprintf('starting slice timing correction run %g subject %g \n',frun,s)
            spm_slice_timing(filesin, sliceorder, refslice, timing,'st_');
        end
        
        % ----------------------
        % Field Map - compute VDM
        % ----------------------
        
        if ~isempty(BIDS.subjects(s).fmap)
            fprintf('computing voxel displacement map %g subject %g \n',frun,s)
            if strcmp(options.despike,'on')
                avg = [filepath filesep 'spmup_mean_st_despiked_' filename ext];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(avg,'file'))
                    spmup_bascis([filepath filesep 'despiked_st_' filename ext],'mean'); % use the mean despiked slice timed EPI for QC
                end
            else
                avg = [filepath filesep 'spmup_mean_st_' filename ext];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(avg,'file'))
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
            
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(vdm{frun},'file'))
                % input images
                if isfield(BIDS.subjects(s).fmap{frun},'phasediff')
                    % two magnitudes (use only 1) and 1 phase difference image
                    matlabbatch = get_FM_workflow('phasediff');
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = {[BIDS_dir filesep subjects{s}.fieldmap(frun).phasediff]};
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {[BIDS_dir filesep subjects{s}.fieldmap(frun).mag1]};
                else
                    % two magnitide images and 2 phase images
                    matlabbatch = get_FM_workflow('phase&mag');
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase = {[BIDS_dir filesep subjects{s}.fieldmap.phase1]};
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag   = {[BIDS_dir filesep subjects{s}.fieldmap.mag1]};
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase  = {[BIDS_dir filesep subjects{s}.fieldmap.phase2]};
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag    = {[BIDS_dir filesep subjects{s}.fieldmap.mag2]};
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
                    if strcmp(BIDS.subjects(s).fmap{frun}.meta.PhaseEncodingDirection,'j') || strcmp(BIDS.subjects(s).fmap{frun}.meta.PhaseEncodingDirection,'y')
                        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = 1;
                    elseif strcmp(BIDS.subjects(s).fmap{frun}.meta.PhaseEncodingDirection,'-j') || strcmp(BIDS.subjects(s).fmap{frun}.meta.PhaseEncodingDirection,'-y')
                        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = -1;
                    end
                elseif isfield(BIDS.subjects(s).func(frun).meta,'PhaseEncodingDirection')
                    if strcmp(BIDS.subjects(s).func(frun).meta.PhaseEncodingDirection,'j') || strcmp(BIDS.subjects(s).func(frun).meta.PhaseEncodingDirection,'y')
                        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = 1;
                    elseif strcmp(BIDS.subjects(s).func(frun).meta.PhaseEncodingDirection,'-j') || strcmp(BIDS.subjects(s).func(frun).meta.PhaseEncodingDirection,'-y')
                        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = -1;
                    end
                else
                    error('No phase encoding direction found')
                end
                
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {avg}; % use the mean despiked / slice timed image
                spm_jobman('run',matlabbatch); clear matlabbatch;
            end
        end
    end
    
    % ------------
    % Realignment
    % ------------
    if isempty(BIDS.subjects(s).fmap)
        % file out of realign will be
        [filepath,filename,ext] = fileparts(st_files{1});
        mean_realigned_file = [filepath filesep 'mean' filename ext];
        for frun = 1:size(BIDS.subjects(s).func,2)
            realigned_file{frun} = st_files{frun}; % because we don't reslice, simple encode the linear transform in the header
            [filepath,filename,ext] = fileparts(st_files{frun});
            multi_reg{frun} = [filepath filesep 'rp_' filename '.txt'];
        end
        
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(mean_realigned_file,'file'))
            fprintf('subject %g: starting realignment \n',s)
            for frun = 1:size(BIDS.subjects(s).func,2)
                matlabbatch{1}.spm.spatial.realign.estwrite.data{frun} = {st_files{frun}};
            end
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp  = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight  = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which   = [0 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp  = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask    = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';
            spm_jobman('run',matlabbatch); clear matlabbatch;
            disp('setting up AC coordinate based on template')
            tmp = realigned_file; tmp{length(tmp)+1} = mean_realigned_file;
            spmup_auto_reorient(tmp',length(tmp)); clear tmp;
        end
    else
        fprintf('subject %g: starting realignment and unwarping \n',s)
        % file out of realign will be
        [~,filename,ext] = fileparts(st_files{1});
        mean_realigned_file = [filepath filesep 'meanur' filename ext];
        for frun = 1:size(BIDS.subjects(s).func,2)
            [~,filename,ext] = fileparts(st_files{frun});
            realigned_file{frun} = [filepath filesep 'ur' filename ext]; % because we have the reslice here (not linear)
            multi_reg{frun} = [filepath filesep 'rp_' filename '.txt'];
        end
        
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(mean_realigned_file,'file'))
            for frun = 1:size(BIDS.subjects(s).func,2)
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
            disp('setting up AC coordinate based on template')
            tmp = realigned_file; tmp{length(tmp)+1} = mean_realigned_file;
            spmup_auto_reorient(tmp',length(tmp)); clear tmp
        end
    end
    
    % -----------------------------------------------
    % Coregistration / Segmentation / Normalization
    % ----------------------------------------------
    [filepath,filename,ext] = fileparts(subjects{s}.anat);
    normalization_file = [BIDS.dir filesep filepath filesep 'y_' filename ext];
    normalizated_file  = [BIDS.dir filesep filepath filesep 'wm' filename ext];
    class{1} = [BIDS.dir filesep filepath filesep 'c1' filename ext];
    class{2} = [BIDS.dir filesep filepath filesep 'c2' filename ext];
    class{3} = [BIDS.dir filesep filepath filesep 'c3' filename ext];
    EPI_class{1} = [BIDS.dir filesep filepath filesep 'c1r' filename ext];
    EPI_class{2} = [BIDS.dir filesep filepath filesep 'c2r' filename ext];
    EPI_class{3} = [BIDS.dir filesep filepath filesep 'c3r' filename ext];
    Normalized_class{1} = [BIDS.dir filesep filepath filesep 'wc1' filename ext];
    Normalized_class{2} = [BIDS.dir filesep filepath filesep 'wc2' filename ext];
    Normalized_class{3} = [BIDS.dir filesep filepath filesep 'wc3' filename ext];
    for frun = 1:size(BIDS.subjects(s).func,2)
        [filepath,filename,ext] = fileparts(realigned_file{frun});
        Normalized_files{frun} = [filepath filesep 'w' filename ext];
    end
    
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(normalization_file,'file'))
        fprintf('subject %g: coregister, segment and normalize \n',s)
        if exist('matlabbatch','var')
            clear matlabbatch
        end
        % coregister anotomical to mean EPI
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {mean_realigned_file};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {subjects{s}.anat};
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        % segment the coregistered T1
        matlabbatch{2}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate & Reslice: Coregistered Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
        matlabbatch{2}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{2}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{2}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{2}.spm.spatial.preproc.tissue(1).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,1']};
        matlabbatch{2}.spm.spatial.preproc.tissue(1).ngaus = 2;
        matlabbatch{2}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(2).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,2']};
        matlabbatch{2}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{2}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(3).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,3']};
        matlabbatch{2}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{2}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(4).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,4']};
        matlabbatch{2}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{2}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(5).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,5']};
        matlabbatch{2}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{2}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(6).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,6']};
        matlabbatch{2}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{2}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{2}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{2}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{2}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{2}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{2}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{2}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{2}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{2}.spm.spatial.preproc.warp.write = [1 1];
        % normalize T1 and EPI data
        matlabbatch{3}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
        matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
        matlabbatch{3}.spm.spatial.normalise.write.subj.resample(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
        matlabbatch{3}.spm.spatial.normalise.write.subj.resample(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
        matlabbatch{3}.spm.spatial.normalise.write.subj.resample(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
        matlabbatch{3}.spm.spatial.normalise.write.subj(2).def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
        for frun = 1:size(BIDS.subjects(s).func,2)
            matlabbatch{3}.spm.spatial.normalise.write.subj(2).resample(frun,:) = {realigned_file{frun}};
        end
        matlabbatch{3}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70 ; 78 76 85];
        matlabbatch{3}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{3}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{3}.spm.spatial.normalise.write.woptions.prefix = 'w';
        % re-run the segmentation on the resliced T1 to get tissue classes in
        % the same space as the EPI data before mormalization
        matlabbatch{4}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
        matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{4}.spm.spatial.preproc.channel.write = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,1']};
        matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus = 2;
        matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,2']};
        matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,3']};
        matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,4']};
        matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,5']};
        matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,6']};
        matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{4}.spm.spatial.preproc.warp.write = [0 0];
        spm_jobman('run',matlabbatch); clear matlabbatch;
    end
    
    % --------------
    % Smoothing
    % --------------
    if strcmp(options.derivatives,'off') % otherwise do it after stats
        for frun = 1:size(BIDS.subjects(s).func,2)
            [filepath,filename,ext] = fileparts(Normalized_files{frun});
            stats_ready{frun} = [filepath filesep 's' filename ext];
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(stats_ready{frun},'file'))
                fprintf('subject %g: smoothing run %g \n',s,frun);
                spm_smooth(Normalized_files{frun},stats_ready{frun},options.skernel);
            end
        end
    else
        stats_ready = Normalized_files;
    end
    
    % ----------------------------
    % QC and additional regressors
    % ----------------------------
    if strcmp(options.QC,'on') %
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist([fileparts(normalizated_file) filesep 'anatQA.mat'],'file'))
            fprintf('subject %g: Anatomical Quality control \n',s)
            % Basic QA for anatomical data is to get SNR, CNR, FBER and Entropy
            % This is useful to check coregistration and normalization worked fine
            % Doing the analysis pre- post- processing also allows to see how much
            % deformation / changes were introduced and if we have groups this is
            % definitively something the check
            anatQA{s}.subject_space  = spmup_anatQA([BIDS.dir filesep subjects{s}.anat],class{1},class{2});
            anatQA{s}.template_space = spmup_anatQA(normalizated_file,Normalized_class{1},Normalized_class{2});
            tmp = anatQA{s}; save([fileparts(normalizated_file) filesep 'anatQA.mat'],'tmp'); clear tmp
        end
        
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist([fileparts(normalizated_file) filesep 'anatQA.mat'],'file'))
            fprintf('subject %g: fMRI Quality control \n',s)
            % For functional data, QA is consists in getting temporal SNR and then
            % check for motion - here we also compute additional regressors to
            % account for motion
            if s==1 && exist([pwd filesep 'tSNR.ps'],'file')
                delete([pwd filesep 'tSNR.ps'])
            end
            
            for frun = 1:size(BIDS.subjects(s).func,2)
                figure('Name','title')
                set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
                title(['subject ' num2str(s) ' run ' num2str(frun)],'FontSize',20);  box off;
                set(gca,'xcolor',get(gcf,'color')); set(gca,'xtick',[]);
                set(gca,'ycolor',get(gcf,'color')); set(gca,'ytick',[]);drawnow;
                if exist([pwd filesep 'tSNR.ps'],'file')
                    print (gcf,'-dpsc2', '-bestfit', [pwd filesep 'tSNR.ps'], '-append');
                else
                    print (gcf,'-dpsc2', '-bestfit', [pwd filesep 'tSNR.ps']);
                end
                close('title')
                tSNR{s,frun} = spmup_temporalSNR(stats_ready{frun},EPI_class,1);
                tmp = tSNR{s}; save([fileparts(normalizated_file) filesep 'tSNR.mat'],'tmp'); clear tmp
            end
        end
    end
    
    
    % --------------------
    % 1st level modelling
    % --------------------
    filepath = fileparts(fileparts([BIDS_dir filesep subjects{s}.func(frun,:)]));
    SPMmat_file = [filepath filesep 'Stats' filesep 'SPM.mat'];
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(SPMmat_file,'file'))
        
        mkdir([filepath filesep 'Stats'])
        matlabbatch{1}.spm.stats.fmri_spec.dir = {[filepath filesep 'Stats']};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = BIDS.subjects(1).func(1).meta.RepetitionTime;
        N = length(BIDS.subjects(1).func(1).meta.SliceTiming);
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = N;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(N/ 2);
        
        % sessions, onsets and durations
        for frun = 1:size(BIDS.subjects(s).func,2)
            N_events = length(BIDS.subjects(s).func(frun).events.trial_type);
            for n=1:N_events
                cond{n} = BIDS.subjects(s).func(frun).events.trial_type(n,:);
            end
            cond = unique(cond);
            N_cond = length(cond);
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans = {stats_ready{frun}};
            onsets = NaN(N_events,N_cond);
            durations = NaN(N_events,N_cond);
            for C = 1:N_cond
                for n=1:N_events
                    if BIDS.subjects(s).func(frun).events.trial_type(n,:) == cond{C}
                        onsets(n,C) = BIDS.subjects(s).func(frun).events.onset(n);
                        durations(n,C) = BIDS.subjects(s).func(frun).events.duration(n);
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
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        save([filepath filesep 'stats_batch_sub' num2str(s) '.mat'],'matlabbatch');
        spm_jobman('run',matlabbatch); clear matlabbatch;
    end
end

disp('spmup BIDS done - all subjects processed')

end

%% helper functions
% -------------------------------------------------------------------------
function options = get_all_options(BIDS_dir)
% routine that returns all available options
% currently set to run the minimal pipeline

options.keep_data ='off';
options.overwrite_data ='on';
options.removeNvol = 0;
options.outdir = [BIDS_dir filesep 'spmup_BIDS_processed'];

% depiking
options.despike= 'off';
options.despiking_window = [];
% realign
options.ignore_fieldmaps = 'on';
% Drifter
options.drifter= 'off';
% Noise regressors
options.motionexp = 'off'; % 24 motion param
options.scrubbing = 'off'; % scrubbing
options.compcor= 'off'; % compcor
% Smoothing
options.skernel = [6 6 6];
% 1st level analysis
options.stats = 'on';
options.derivatives = 'off';
% QC
options.QC = 'on';

end

% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function job = get_first_level_param


end
