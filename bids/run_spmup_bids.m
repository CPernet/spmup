function [subjects,opt] = run_spmup_bids(BIDS, subjects, varargin)

% SPMUP can be run using SPM batch mode but you can also analyze an entire
% fMRI BIDS dataset automatically using this function - it will unpack zip
% files, preprocess (despike, slice time, realign, smooth), QA (anat, fMRI),
% save the augmented regressors for motion, censoring, etc, and run subject
% level GLM (augmented regressors, amplitude correction, delay estimation )
% using almost all the gizmos of the toolbox. Note the GLM will be the
% simplest possible version, i.e. 1 regressor per unique events, which
% might not be suitable, although using contrasts any other combinations
% can, in theory, be obtained.
%
% FORMAT run_spmup_bids(BIDS, subjects, options)
%
% INPUT BIDS the dataset properties
%       subjects is the structure coming out of spmup_BIDS_unpack
%       options is either a set of key/value pairs matching fields from
%               spmup_getoptions or direction the options structure
%
% OUTPUT the preprocessed and analyzed data are written on disk
%        subjects is the updated structure with processed file names, QC, etc
%        opt is a strcuture of options per subject
%        (subjects and opt are also saved on disk)
%
% Example usage       BIDS_dir        = '/data1/cpernet/ds000117';
%                     options         = spmup_getoptions(BIDS_dir);
%                     options.anat    = {'T1w'};
%                     [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options)
% you may need to run system(['chmod -Rf 755 ' options.outdir]) on linux servers
%                     subjects = run_spmup_bids(BIDS,subjects,options)
%
% Cyril Pernet July 2022
% --------------------------
%  Copyright (C) SPMUP Team

%% Set options, matlab path and get data

% make sure all folders are in the path
addpath(genpath(fullfile(fileparts(which('spm.m')),['toolbox' filesep 'spmup'])));
if nargin == 0
    help run_spmup_bids
    return
end

% check inputs
if nargin >= 3
    if nargin == 3 && isstruct(varargin{1})
        options = varargin{1};
    else
        for v=1:2:length(varargin)
            if any(contains(varargin{v},fieldnames(options)))
                options.(varargin{v}) = varargin{v+1};
            end
        end
    end
end

% get the initial subject BIDS structure and check task(s)
subjs_ls  = spm_BIDS(BIDS,'subjects');
if isempty(options.task)
    % assume all task the same
    clear task
    index = 1;
    for s=1:size(BIDS.subjects,2)
        if ~isempty(BIDS.subjects(s).func)
            for r=1:size(BIDS.subjects(s).func,2)
                task{index} = BIDS.subjects(s).func(r).task; %#ok<AGROW>
                index = index +1;
            end
        end
    end
    if all(size(unique(task))==[1 1])
        options.task = {task{1}};
        Ntask        = 1;
        clear task;
    else
        error('multiple tasks detected, please specify with options what to do')
    end
else
    Ntask = length(options.task);
end

% replicate options per subject
opt = get_subject_options(options,subjects);

%% process

spmup_setparallel(options.Ncores)

for task = 1:Ntask
    for s = 1:numel(subjects)
        if ~isempty(options.task)
            opt(s).task = options.task{task};
        end
        sess_func   = size(subjects{s}.func,1);
        sess_anat   = size(subjects{s}.anat,1);
        sess_names  = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
        if isempty(sess_names)
            sess_names = {''};
        end
        
        for session = 1:sess_func
            
            if ~isempty(options.task)
                run_ls = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
                    'ses', sess_names{session}, 'task', opt(s).task, 'type', 'bold');
            else
                run_ls = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
                    'ses', sess_names{session}, 'type', 'bold');
            end
            
            if ~isempty(run_ls)
                metadata = cell(size(run_ls));
                % in the new structure, we added run folder to process data
                % this needs to be matched as BIDS var does not have this
                for r=1:size(run_ls,1)
                    
                    % anat
                    if sess_anat == 1 && session > 1 % only one anat for multiple sesions
                        subjects{s}.anat{session,:} = subjects{s}.anat{1,:}; % replicate
                    end
                    subject_sess.anat               = subjects{s}.anat(session,cellfun(@(x) ~isempty(x), subjects{s}.anat(session,:)))';
                    if ~iscell(subject_sess.anat)
                        subject_sess.anat           = {subject_sess.anat};
                    end
                    
                    % fmri
                    % BIDS_name: % e.g. 'sub-53888_ses-02_acq-ep2d_task-rest_bold'
                    [~,BIDS_name, ext] = fileparts(run_ls{r}); 
                    if strcmp(ext, '.gz')
                        [~,BIDS_name] = fileparts(BIDS_name);
                    end
                    % sub_names : e.g. {'sub-53888_ses-02_acq-ep2d_task-faces_bold'}    {'sub-53888_ses-02_acq-ep2d_task-rest_bold'}
                    tmp = subjects{s}.func(session,:)';
                    for f = numel(tmp):-1:1
                        [~,sub_names{f,1}] = fileparts(tmp{f});
                    end
                    
                    % match BIDS_name with the right run in the subjects structure
                    run_index     = find(strcmpi(sub_names,BIDS_name)); 
                    if isempty(run_index) % because shit happens
                        try
                            for subr = 1:size(sub_names,2)
                                run_index(subr) = strncmp(sub_names{subr},BIDS_name(1:length(sub_names{subr})),length(sub_names{subr}));
                            end
                            run_index = find(run_index);
                        catch
                            run_index = 1; warning('multiple runs, but indexing failed - using simple ordering')
                        end
                    end
                    
                    if isempty(run_index)
                        error('somehow we cannot match BIDS_name with the spmup unpacked subject structure information, some silly matching error to figure out - please reach out/raise an issue so we can fix it')
                    else
                        run_ls{r}     = cell2mat(subjects{s}.func(session,run_index)); % which gives the right path
                        metadata{r}   = cell2mat(subjects{s}.func_metadata(session,run_index)); % and thus metadata
                    end
                end
                
                subject_sess.func          = run_ls;
                subject_sess.func_metadata = metadata;
                subject_sess.fieldmap      = {};
                if isfield(subjects{s}, 'fieldmap')
                    subject_sess.fieldmap  = subjects{s}.fieldmap(session);
                end
                if isfield(subjects{s},'event')
                    event_index            = contains(subjects{s}.event(session,:),opt(s).task);
                    subject_sess.event     = subjects{s}.event(session,event_index)';
                end
                
                % compute
                [subject_sess,opt(s)]      = spmup_BIDS_preprocess(subject_sess, opt(s));
                if strcmpi(opt(s).GLM,'on')
                    subject_sess = spmup_BIDS_1rstlevel(subject_sess, opt(s));
                end
                
                % update subjects structure
                subjects{s}.anat{session,1} = subject_sess.anat{1};
                for run = 1:size(subjects{s}.func,2)
                    [~,tmp] = fileparts(subjects{s}.func{session,run});
                    if contains(tmp,opt(s).task)
                        current_run_task = tmp(strfind(tmp,'run-'):strfind(tmp,'_bold')-1);
                        for srun = 1:max(size(subject_sess.func))
                            if contains(subject_sess.func{srun},current_run_task)
                                subjects{s}.func{session,run} = subject_sess.func{srun};
                            end
                        end
                    end
                end
                clear run
                subjects{s}.motionfile{session} = subject_sess.motionfile;
                subjects{s}.tissues{session}    = subject_sess.tissues;
                subjects{s}.func_qa{session}    = subject_sess.func_qa;
                subjects{s}.anat_qa{session}    = subject_sess.anat_qa;
                if isfield(subject_sess,'stats')
                    subjects{s}.stats{session}  = subject_sess.stats;
                end
                if isfield(subject_sess,'boostedcontrast')
                    subjects{s}.boostedcontrast = subject_sess.boostedcontrast;
                end
                
                clear run_index event_index subject_sess
                save([options.outdir filesep 'spmup_subjects_task-' options.task{task} '.mat'],'subjects');
                save([options.outdir filesep 'spmup_options_task-'  options.task{task} '.mat'],'opt');
            end
        end
        
        disp('---------------------------------------')
        fprintf('session %s finished!\n',sess_names{session})
        disp('---------------------------------------')
        
        %% reinitialize opt for the next session
        opt = get_subject_options(options,subjects);
    end
    
    %% save and report QC
    if iscell(sess_names)
        has_session = ~isempty(sess_names{1});
    else
        has_session = ~isempty(sess_names);
    end
    
    table_name  = spmup_BIDS_QCtables(subjects, 'anat');
    for t=1:size(table_name,2)
        basename = spm_file(table_name{t}, 'basename');
        if ~has_session
            destination = 'AnatQC';
        else
            destination = basename;
        end
        destination = fullfile(options.outdir, [destination, '_task-' options.task{task}, '.tsv']);
        % fprintf('\nMoving QC results\n\tfrom:%s\n\tto:%s\n', table_name{t}, destination);
        if exist(table_name{t}, 'file') && isfolder(fileparts(destination))
            movefile(table_name{t},destination);
        end
    end
    
    table_name  = spmup_BIDS_QCtables(subjects, 'fMRI');
    for t=1:size(table_name,2)
        basename = spm_file(table_name{t}, 'basename');
        has_run = ~isempty(strfind(basename,'run'));
        if ~has_session && ~has_run
            destination = 'fMRIQC';
        else
            destination = basename;
        end
        destination = fullfile(options.outdir, [destination '_task-' options.task{task} '.tsv']);
        % fprintf('\nMoving QC results\n\tfrom:%s\n\tto:%s', table_name{t}, destination);
        if exist(table_name{t}, 'file') && isfolder(fileparts(destination))
            movefile(table_name{t},destination);
        end
    end
    
    disp('---------------------------------------')
    fprintf('task %s finished!\n',options.task{task})
    disp('---------------------------------------')
end


function opt = get_subject_options(options,subjects)

% routine to replicate options per subject
if ~isempty(options.VDM)
    if numel(subjects) ~= length(options.VDM)
        error('the number of subjects and the VDM cell array in options must be the same size (if multiple VDMs per subject, use a cell array of cells)')
    end
    for s=length(options.VDM):-1:1
        opt(s) = options;
        opt(s).VDM = options.VDM{s}; % allows any number of FM/VDM
    end
else
    for s = numel(subjects):-1:1
        opt(s) = options; % .VDM needs updated if empty anyway
    end
end
