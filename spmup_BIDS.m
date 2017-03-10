function spmup_BIDS(BIDS_dir,options)

% routine to read and unpack BIDS fMRI
% does a fixed data preprocessing and 1st level analysis
% 1.  depiking
% 2.  slice timing
% 3.  AC coordinate setting
% 4a. Computing VDM maps (if field maps)
% 4b. Realignment (with unwarping if 4a)
% 5.  Drifter (removes cardiac and respiratory signals)
% 6.  Coregister T1 to EPI, segment and normalize it
% 7.  Normalize EPI data
% 8.  Generate QA report and noise regressors
% 9.  Performs 1st level analysis with 1st derivative
% 10. Boost beta estimates and smooth them
%
% options.outdir = where to write the data
% options.preprocess.removeNvol = number of initial volumes to remove
% options.keep_data 'off' (default) or 'on' to keep all steps - off means
%                    only the last processed data are available (typically
%                    despiked, slice timed, realigned, normalized data)
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

if nargin == 0 || nargin == 1
    options = get_all_options(BIDS_dir);
end

if ~isfield(options.preprocess,'removeNvol') || isempty(options.preprocess.removeNvol)
    options.preprocess.removeNvol = input('How many initial EPI volumes to discard?: '); disp(' ')
end

spm('defaults', 'FMRI');

%% get the data info
% -------------------------------------------------------------------------
cd(BIDS_dir)
BIDS = spm_BIDS(BIDS_dir);

%% unpack data
% -------------------------------------------------------------------------
if ~isfield(options,'outdir')
    options.outdir = [BIDS_dir filesep 'spmup_BIDS_processed'];
end

mkdir(options.outdir);
for s=1:size(BIDS.subjects,2)
    fprintf('subject %g: unpacking anatomical data \n',s)
    % when encountering longitudinal dataset this needs updated for 'session'
    in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat.filename];
    out = [options.outdir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat.filename(1:end-3)];
    gunzip(in, out); spmup_auto_reorient([out filesep spm_file(out,'basename') '.nii']); % here set AC coordinate
end

disp('unpacking functional data')
if size(BIDS.subjects(s).func,2) >= 2
    parpool('local',feature('numCores')-1); % use all available cores -1
end

for s=1:size(BIDS.subjects,2)
    fprintf('subject %g: unpacking functional data \n',s)
    parfor frun = 1:size(BIDS.subjects(s).func,2)
        in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'func' filesep BIDS.subjects(s).func(frun).filename];
        out = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep BIDS.subjects(s).func(frun).filename(1:end-3)];
        gunzip(in, out);
    end
end

if ~isempty(BIDS.subjects(1).fmap)
    disp('unpacking field maps')
end

parfor s=1:size(BIDS.subjects,2)
    % this needs update for longitinal dataset ie multiple 'sessions'
    for frun = 1:size(BIDS.subjects(s).fmap,2)
        for f=1:size(BIDS.subjects(s).fmap(frun).filenames,1)
            in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap(frun).filenames{f}];
            out = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap(frun).filenames{f}(1:end-3)];
            gunzip(in, out);
        end
    end
end

try delete(gcp('nocreate')); end

%% run preprocessing using options
% -------------------------------------------------------------------------
spm_jobman('initcfg');

for s=1:size(BIDS.subjects,2)
    fprintf('subject %g: starting despiking and slice timing correction for each run\n',s)
    for frun = 1:size(BIDS.subjects(s).func,2)
        % deal with uncompressed data
        [~,file,~]=fileparts(BIDS.subjects(s).func(frun).filename);
        pth = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep file];
        filesin = [pth filesep file];
        despiked_files = [pth filesep 'despiked_' file];
        st_files{frun} = [pth filesep 'st_despiked_' file];
        avg = [pth filesep 'spmup_mean_st_despiked_' file];
        
        % -----------------------------------------
        % despiking using adaptive median filter
        % -----------------------------------------
        flags = struct('auto_mask','on', 'method','median', 'window', ...
            options.preprocess.despiking.window,'skip',options.preprocess.removeNvol);
        cd(pth); spmup_despike(filesin,[],flags);
        if strcmpi(options.keep_data,'off')
            delete(filesin);
        end
        cd(BIDS_dir);
        
        % -------------
        % AC coordinate
        % -------------
        spmup_auto_reorient(despiked_files,1);
        
        % -------------
        % slice timing
        % -------------
        % sliceorder - vector containig the acquisition time for each slice in milliseconds
        % refslice   - time in milliseconds for the reference slice
        % timing     - [0 TR]
        
        sliceorder  = BIDS.subjects(s).func(frun).meta.SliceTiming; % time 
        refslice    = sliceorder(round(length(BIDS.subjects(s).func(frun).meta.SliceTiming)/2)); % time
        timing      = [0 BIDS.subjects(s).func(frun).meta.RepetitionTime];
        spm_slice_timing(despiked_files, sliceorder, refslice, timing,'st_');
        
        % ----------------------
        % Field Map - compute VDM
        % ----------------------
        
        if ~isempty(BIDS.subjects(s).fmap(frun).filenames)
            disp('computing voxel displacement map ')
            spmup_bascis(st_files{frun},'mean'); % use the mean EPI (before realign as ref)
            
            % input images
            if size(BIDS.subjects(s).fmap(frun).filenames,1) <= 3
                % two magnitudes and 1 phase difference image
                matlabbatch = get_FM_workflow('phasediff');
                for ii=1:size(BIDS.subjects(s).fmap(frun).filenames,1)
                    phase_img(ii) = ~isempty(findstr(BIDS.subjects(s).fmap(frun).filenames{ii},'phase'));
                end
                [~,file,~]=fileparts(BIDS.subjects(s).fmap(frun).filenames{find(phase_img)});
                pth = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps'];
                phase_img = [pth filesep file filesep file];
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = {phase_img};
                magnitude_img = [pth filesep BIDS.subjects(s).fmap(frun).filenames{1}(1:end-3) filesep BIDS.subjects(s).fmap(frun).filenames{1}(1:end-3)];
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {magnitude_img};
            else
                % two magnitide images and 2 phase images
                matlabbatch = get_FM_workflow('phase&mag');
                for ii=1:size(BIDS.subjects(s).fmap(frun).filenames,1)
                    phase_img(ii) = ~isempty(findstr(BIDS.subjects(s).fmap(frun).filenames{ii},'phase'));
                end
                pth = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps'];
                phase_index = find(phase_img);
                [~,file,~]=fileparts(BIDS.subjects(s).fmap(frun).filenames{phase_index(1)});
                phase_img = [pth filesep file filesep file];
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase = {phase_img};
                magnitude_img = [pth filesep BIDS.subjects(s).fmap(frun).filenames{1}(1:end-3) filesep BIDS.subjects(s).fmap(frun).filenames{1}(1:end-3)];
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag = {magnitude_img};
                [~,file,~]=fileparts(BIDS.subjects(s).fmap(frun).filenames{phase_index(2)});
                phase_img = [pth filesep file filesep file];
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase = {phase_img};
                magnitude_img = [pth filesep BIDS.subjects(s).fmap(frun).filenames{2}(1:end-3) filesep BIDS.subjects(s).fmap(frun).filenames{2}(1:end-3)];
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag = {magnitude_img};
            end
            
            % update parameters
            echotimes =  1000.*[BIDS.subjects(s).fmap(frun).meta.EchoTime1 BIDS.subjects(s).fmap(frun).meta.EchoTime2]; % change to ms
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = echotimes;
            
            if isfield(BIDS.subjects(s).fmap(frun).meta,'TotalReadoutTime')
                TotalReadoutTime = BIDS.subjects(s).fmap(frun).meta.TotalReadoutTime;
            elseif isfield(BIDS.subjects(s).func(frun).meta,'RepetitionTime')
                TotalReadoutTime = BIDS.subjects(s).fmap(frun).meta.RepetitionTime;
            elseif isfield(BIDS.subjects(s).func(frun).meta,'EffectiveEchoSpacing')
                TotalReadoutTime = (BIDS.subjects(s).func(frun).meta.NumberOfEchos-1)*BIDS.subjects(s).func(frun).meta.EffectiveEchoSpacing;
            end
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = TotalReadoutTime;
            
            if isfield(BIDS.subjects(s).fmap(frun).meta,'PulseSequenceType')
                if sum(findstr(BIDS.subjects(s).fmap(frun).meta.PulseSequenceType,'EPI')) ~= 0
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
                else
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
                end
            else
                 disp('using default sequence! assuming non-EPI acquisition')
            end
            
            if isfield(BIDS.subjects(s).fmap(frun).meta,'PhaseEncodingDirection')
                if strcmp(BIDS.subjects(s).fmap(frun).meta.PhaseEncodingDirection,'j') || strcmp(BIDS.subjects(s).fmap(frun).meta.PhaseEncodingDirection,'y')
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = 1;
                elseif strcmp(BIDS.subjects(s).fmap(frun).meta.PhaseEncodingDirection,'-j') || strcmp(BIDS.subjects(s).fmap(frun).meta.PhaseEncodingDirection,'-y')
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

            
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {avg}; % use the mean despiked / AC oriented / slice timed image
            clear phase_img ; out = spm_jobman('run',matlabbatch);
            vdm{frun} = cell2mat(out{1}.vdmfile{1}); cd(BIDS_dir);
        end
    end
    
    % ------------
    % Realignment
    % ------------
    fprintf('subject %g: starting realignment and unwarping (if field maps available)\n',s)
    
    if isempty(BIDS.subjects(s).fmap(1).filenames)
        matlabbatch{1}.spm.spatial.realign.estwrite.data = st_files;
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
    else
        for frun = 1:size(BIDS.subjects(s).func,2)
            matlabbatch{1}.spm.spatial.realignunwarp.data(frun).scans = st_files{frun};
            matlabbatch{1}.spm.spatial.realignunwarp.data(frun).pmscan = vdm{frun};
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
    end
    R{s} = spm_jobman('run',matlabbatch); 
     
    
    
    % ----------------------------
    % QA and additional regressors
    % ----------------------------
    
    
    
    % -----------------------
    % Coregistration with T1
    % ----------------------
    
    
    % --------------
    % Normalization
    % --------------
    
    
    % --------------
    % Smoothing
    % --------------
    
    
end


% --------------------
% 1st level modelling
% --------------------



disp('spmup BIDS done - all subjects processed')

end

%% helper functions
% -------------------------------------------------------------------------
function options = get_all_options(BIDS_dir)
% routine that returns all available options

preprocess_options.removeNvol = [];
% despiking
preprocess_options.despiking.window = [];
preprocess_options.drifter = [];

options = struct('outdir',[BIDS_dir filesep 'spmup_BIDS_processed'], ...
    'preprocess',preprocess_options,'keep_data','off');
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
function job = get_realignment_param


end

% -------------------------------------------------------------------------
function job = get_first_level_param


end
