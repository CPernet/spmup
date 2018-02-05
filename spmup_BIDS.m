function [anatQA,fMRIQA]=spmup_BIDS(BIDS_dir,choices)

% routine to read and unpack BIDS fMRI and preprocess data, build the first
% level analysis - with various options available
%
% FORMAT spmup_BIDS
%        spmup_BIDS(BIDS_dir)
%        spmup_BIDS(BIDS_dir,choices)
%
% INPUTS BIDS_dir is the BIDS directory
%        choices a structure with the following fields:
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
%     'ignore_fieldmaps', 'on',  'outdir', 'spmup_BIDS_processed_basic', 'QC', 'off'); % standard SPM pipeline
% [anatQA,fMRIQA]=spmup_BIDS(pwd,choice)
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

% inputs
% -------------------------------------------------------------------------
if nargin == 0
    BIDS_dir= uigetdir(pwd,'select BIDS directory');
    if BIDS_dir == 0
        return
    else
        choices = get_all_options(BIDS_dir);
    end
end

% get default options
% --------------------
options = get_all_options(BIDS_dir);

% update options with user input
% --------------------------------
I=intersect(fieldnames(options),fieldnames(choices));
for f=1:length(I)
    options = setfield(options,I{f},getfield(choices,I{f}));
end
clear choices
spm('defaults', 'FMRI');

% get the data info
% --------------------
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

if ~exist(options.outdir,'dir')
    mkdir(options.outdir);
end

if isempty(BIDS.sessions)
    BIDS.sessions = 1;
end

% unpack anat and center [0 0 0]
% -------------------------------
for s=1:size(BIDS.subjects,2) % for each subject
    for session = 1:BIDS.sessions % for each session
        in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat(session).filename];
        if strcmp(in(end-2:end),'.gz') % if compressed
            subjects{s}.anat = [options.outdir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat(session).filename(1:end-3)];
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.anat,'file'))
                fprintf('subject %g: unpacking anatomical data \n',s)
                gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'anat' ]);
                spmup_auto_reorient(subjects{s}.anat); disp('anat reoriented');
            end
        else % if nor compressed
            subjects{s}.anat = [options.outdir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat(session).filename];
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.anat,'file'))
                fprintf('subject %g: copying anatomical data \n',s)
                copyfile(in, [options.outdir filesep BIDS.subjects(s).name filesep 'anat' ]);
                spmup_auto_reorient(subjects{s}.anat); disp('anat reoriented');
            end
        end
    end
end

% if mote than 2 functional run, try to unpack data using multiple cores
% ----------------------------------------------------------------------
if size(BIDS.subjects(s).func,2) >= 2
    try
        parpool(feature('numCores')-1); % use all available cores -1
    catch no_parpool
        disp(no_parpool.message)
    end
end

% unpack functional and field maps (still need to work out sessions here)
% -----------------------------------------------------------------------
parfor s=1:size(BIDS.subjects,2)
    
    for frun = 1:size(BIDS.subjects(s).func,2) 
        
        % functional
        in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'func' filesep BIDS.subjects(s).func(frun).filename];
        if strcmp(in(end-2:end),'.gz')
            subjects{s}.func(frun,:) = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep BIDS.subjects(s).func(frun).filename(1:end-3)];
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.func(frun,:),'file'))
                fprintf('subject %g: unpacking functional data run %g \n',s,frun)
                gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun)]);
            end
        elseif strcmp(in(end-2:end),'nii')
            subjects{s}.func(frun,:) = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep BIDS.subjects(s).func(frun).filename];
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.func(frun,:),'file'))
                fprintf('subject %g: copying functional data \n',s)
                copyfile(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun)]);
            end
        end
        
        % field maps
        if ~isempty(BIDS.subjects(s).fmap)
            fprintf('subject %g: checking field maps \n',s)
            in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap{frun}.magnitude1];
            if strcmp(in(end-2:end),'.gz')
                subjects{s}.fieldmap(frun,:).mag1 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.magnitude1(1:end-3)];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).mag1,'file'))
                    gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                end
            elseif strcmp(in(end-2:end),'nii')
                subjects{s}.fieldmap(frun,:).mag1 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.magnitude1];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).mag1,'file'))
                    copyfile(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                end
            end
            
            in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap{frun}.magnitude2];
            if strcmp(in(end-2:end),'.gz')
                subjects{s}.fieldmap(frun,:).mag2 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.magnitude2(1:end-3)];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).mag2,'file'))
                    gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                end
            elseif strcmp(in(end-2:end),'nii')
                subjects{s}.fieldmap(frun,:).mag2 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.magnitude2];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).mag2,'file'))
                    copyfile(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                end
            end
             
             
            if isfield(BIDS.subjects(s).fmap{frun},'phasediff')
                in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap{frun}.phasediff];
                if strcmp(in(end-2:end),'.gz')
                    subjects{s}.fieldmap(frun,:).phasediff = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.phasediff(1:end-3)];
                    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).phasediff,'file'))
                        gunzip(in,[options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                    end
                elseif strcmp(in(end-2:end),'nii')
                    subjects{s}.fieldmap(frun,:).phasediff = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.phasediff];
                    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).phasediff,'file'))
                        copyfile(in,[options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                    end
                end
                
            else
                in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap{frun}.phase1];
                if strcmp(in(end-2:end),'.gz')
                    subjects{s}.fieldmap(frun,:).phase1 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.phase1(1:end-3)];
                    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).phase1,'file'))
                        gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                    end
                else
                    subjects{s}.fieldmap(frun,:).phase1 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.phase1];
                    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).phase1,'file'))
                        copyfile(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                    end
                end
                
                in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'fmap' filesep BIDS.subjects(s).fmap{frun}.phase2];
                if strcmp(in(end-2:end),'.gz')
                    subjects{s}.fieldmap(frun,:).phase2 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.phase2(1:end-3)];
                    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).phase2,'file'))
                        gunzip(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                    end
                else
                    subjects{s}.fieldmap(frun,:).phase2 = [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps' filesep BIDS.subjects(s).fmap{frun}.phase2];
                    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(subjects{s}.fieldmap(frun,:).phase2,'file'))
                        copyfile(in, [options.outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep 'fieldmaps']);
                    end
                end
            end
        end
    end
end

disp('spmup has finished unpacking data')

%% run preprocessing using options
% -------------------------------------------------------------------------
if isfield(BIDS,'skip') %% need to check BIDS spec 
    start_at = BIDS.skip;
else
    start_at = 1;
end

% do the computation if needed
parfor s=1:size(BIDS.subjects,2)
    spmup_BIDSjob(BIDS_dir,BIDS,subjects,s,options,start_at)
end

disp('spmup BIDS done - all subjects processed')
try delete(gcp('nocreate')); end

end

%% option sub-function
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
% Normalization
options.norm = 'EPInorm';
% Smoothing
options.skernel = [8 8 8];
% 1st level analysis
options.stats = 'on';
options.derivatives = 'off';
% QC
options.QC = 'on';

end


