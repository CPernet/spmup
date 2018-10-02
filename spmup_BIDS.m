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
%     'ignore_fieldmaps', 'on',  'outdir', ['..' filesep 'derivatives' filesep 'spmup_BIDS_processed'], 'QC', 'off'); % standard SPM pipeline
% [anatQA,fMRIQA]=spmup_BIDS(pwd,choice)
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

% TO DO
% implement keep_data ?
% implement BIDS.skip

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
    options.outdir = [BIDS_dir filesep '..' filesep 'derivatives' filesep  ...
        'spmup_BIDS_processed']; % ?? isn't that set by the get_all_options subfunction ??
end

if ~exist(options.outdir,'dir')
    mkdir(options.outdir);
end

subjs_ls = spm_BIDS(BIDS,'subjects');

% unpack anat and center [0 0 0]
% -------------------------------
for s=1:size(subjs_ls,2) % for each subject

    sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});

    for session = 1:size(sess_ls,2) % for each session

        if size(sess_ls,2)==1
            sess_fodler = '';
        else
            sess_fodler = ['ses-' sess_ls{session}];
        end

        in = char(spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
            'ses', sess_ls{session}, 'type', 'T1w'));

        [~,name,ext,~] = spm_fileparts(in);

        if strcmp(ext,'.gz') % if compressed
            subjects{s}.anat = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                sess_fodler, 'anat', name); %#ok<*AGROW>

            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                    && ~exist(subjects{s}.anat,'file'))
                fprintf('subject %g: unpacking anatomical data \n',s)
                gunzip(in, fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                sess_fodler, 'anat'));
                spmup_auto_reorient(subjects{s}.anat); disp('anat reoriented');
            end

        elseif strcmp(ext,'.nii') % if not compressed
            subjects{s}.anat = fullfile(options.outdir, ['sub-' subjs_ls(s)], ...
                sess_fodler, 'anat', [name ext]);

            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                    && ~exist(subjects{s}.anat,'file'))
                fprintf('subject %g: copying anatomical data \n',s)
                copyfile(in, fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                sess_fodler, 'anat'));
                spmup_auto_reorient(subjects{s}.anat); disp('anat reoriented');
            end

        elseif isempty(name) % if no valid T1w was found
            warning('No valid T1w file found for subject %s - session %s', subjs_ls{s}, sess_ls{session})
        end

    end

end

% if more than 2 functional run, try to unpack data using multiple cores
% ----------------------------------------------------------------------
if size(spm_BIDS(BIDS, 'runs', 'sub', subjs_ls{s}),2) >= 2 %#ok<*UNRCH>
    try
        parpool(feature('numCores')-1); % use all available cores -1
    catch no_parpool
        disp(no_parpool.message)
    end
end

% unpack functional and field maps (still need to work out sessions here)
% -----------------------------------------------------------------------
% for s=1:size(subjs_ls,2)
parfor s=1:size(subjs_ls,2)
    
    sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
    
    % initialize bold file list so that new runs can be appended with end+1
    % even if they are from different sessions
    subjects{s}.func = [];
    
    for session = 1:size(sess_ls,2) % for each session
        
        if size(sess_ls,2)==1
            sess_fodler = '';
        else
            sess_fodler = ['ses-' sess_ls{session}];
        end
        
        run_ls = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
            'ses', sess_ls{session}, 'type', 'bold');
        
        for frun = 1:size(run_ls,1)
            
            % functional
            in = run_ls {frun,1};
            
            [~,name,ext,~] = spm_fileparts(in);
            
            if strcmp(ext,'.gz')
                subjects{s}.func(end+1,:) = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                    sess_fodler, 'func', ['run' num2str(frun)], name); %#ok<*AGROW>
                
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                        && ~exist(subjects{s}.func(frun,:),'file'))
                    fprintf('subject %g: unpacking functional data run %g \n',s,frun)
                    gunzip(in, fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                        sess_fodler, 'func', ['run' num2str(frun)]));
                end
                
            elseif strcmp(ext,'nii')
                subjects{s}.func(end+1,:) = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                    sess_fodler, 'func', ['run' num2str(frun)], [name ext]); %#ok<*AGROW>
                
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                        && ~exist(subjects{s}.func(frun,:),'file'))
                    fprintf('subject %g: copying functional data \n',s)
                    copyfile(in, fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                        sess_fodler, 'func', ['run' num2str(frun)]));
                end
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
if isfield(BIDS,'skip') %% need to check BIDS spec (???)
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
options.outdir = [BIDS_dir filesep '..' filesep 'derivatives' filesep  ...
    'spmup_BIDS_processed'];

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


