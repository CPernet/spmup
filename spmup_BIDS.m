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


%% unpack anat and center [0 0 0]
% -------------------------------
for s=1:size(subjs_ls,2) % for each subject
    
    sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
    
    for session = 1:size(sess_ls,2) % for each session
        
        if size(sess_ls,2)==1
            sess_folder = '';
        else
            sess_folder = ['ses-' sess_ls{session}];
        end
        
        in = char(spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
            'ses', sess_ls{session}, 'type', 'T1w'));
        
        [~,name,ext,~] = spm_fileparts(in);
        
        if strcmp(ext,'.gz') % if compressed
            subjects{s}.anat = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                sess_folder, 'anat', name); %#ok<*AGROW>
            
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                    && ~exist(subjects{s}.anat,'file'))
                fprintf('subject %g: unpacking anatomical data \n',s)
                gunzip(in, fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                    sess_folder, 'anat'));
                spmup_auto_reorient(subjects{s}.anat); disp('anat reoriented');
            end
            
        elseif strcmp(ext,'.nii') % if not compressed
            subjects{s}.anat = fullfile(options.outdir, ['sub-' subjs_ls(s)], ...
                sess_folder, 'anat', [name ext]);
            
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                    && ~exist(subjects{s}.anat,'file'))
                fprintf('subject %g: copying anatomical data \n',s)
                copyfile(in, fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                    sess_folder, 'anat'));
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
for s=1:size(subjs_ls,2)
    % parfor s=1:size(subjs_ls,2)
    
    sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
    
    % initialize bold file list so that new runs can be appended with end+1
    % even if they are from different sessions
    subjects{s}.func = cell(1);
    
    for session = 1:size(sess_ls,2) % for each session
        
        if size(sess_ls,2)==1
            sess_folder = '';
        else
            sess_folder = ['ses-' sess_ls{session}];
        end
        
        run_ls = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
            'ses', sess_ls{session}, 'type', 'bold');
        
        %% functional
        for frun = 1:size(run_ls,1) % for each run
            
            in = run_ls{frun,1};
            
            [~,name,ext] = spm_fileparts(in);
            
            if strcmp(ext,'.gz') % if compressed
                % we keep track of where the files are stored
                subjects{s}.func{end+1,:} = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                    sess_folder, 'func', ['run' num2str(frun)], name); %#ok<*AGROW>
                
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                        && ~exist(subjects{s}.func{frun},'file'))
                    fprintf('subject %g: unpacking functional data run %g \n',s,frun)
                    gunzip(in, fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                        sess_folder, 'func', ['run' num2str(frun)]));
                end
                
            elseif strcmp(ext,'nii') % if not compressed
                subjects{s}.func{end+1,:} = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                    sess_folder, 'func', ['run' num2str(frun)], [name ext]); %#ok<*AGROW>
                
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                        && ~exist(subjects{s}.func{frun},'file'))
                    fprintf('subject %g: copying functional data \n',s)
                    copyfile(in, fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                        sess_folder, 'func', ['run' num2str(frun)]));
                end
            end
        end
        
        
        %% field maps
        if ismember('fmap', spm_BIDS(BIDS,'modalities', 'sub', subjs_ls{s}, ...
                'ses', sess_ls{session}))
            
            fprintf('subject %g: checking field maps \n',s)
            
            target_dir = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                sess_folder, 'fieldmaps');
            
            % check which types of fieldmpas we are dealing with
            fmap_type_ls = spm_BIDS(BIDS, 'types', 'sub', subjs_ls{s}, ...
                'ses', sess_ls{session}, 'modality', 'fmap');
            
            % goes through each type in case we have several types in that
            % data set
            for ifmap_type_ls = 1:size(fmap_type_ls,1)
                
                fmap_ls = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
                    'ses', sess_ls{session}, 'modality', 'fmap', ...
                    'type', fmap_type_ls(ifmap_type_ls));
                
                % for all the fieldmaps
                for ifmap = 1:size(fmap_ls,1)
                    
                    in = fmap_ls{ifmap,1};
                    [path,name,ext] = spm_fileparts(in);
                    
                    [ext,name,compressed] = iscompressed(ext,name);
                    
                    switch fmap_type_ls{ifmap_type_ls}
                        
                        case 'phasediff'
                            subjects{s}.fieldmap(ifmap,:).phasediff = ...
                                fullfile(target_dir, [name ext]);
                            file_exists = exist(subjects{s}.fieldmap(ifmap,:).phasediff,'file');
                            
                        case 'phase12'
                            if contains(name,'phase1')
                                subjects{s}.fieldmap(ifmap,:).phase1 = ...
                                    fullfile(target_dir, [name ext]);
                                file_exists = exist(subjects{s}.fieldmap(ifmap,:).phase1,'file');
                            elseif contains(name,'phase2')
                                subjects{s}.fieldmap(ifmap,:).phase2 = ...
                                    fullfile(target_dir, [name ext]);
                                file_exists = exist(subjects{s}.fieldmap(ifmap,:).phase2,'file');
                            end
                            
                        case 'fieldmap'
                            warning('Fieldmap type of fielmaps not implemented')
                            
                        case 'epi'
                            warning('EPI type of fielmaps not implemented')
                            
                        otherwise
                            warning('%s is an unsupported type of fieldmap', fmap_type_ls(ifmap_type_ls))
                    end
                    
                    unzip_or_copy(compressed, options, file_exists, in, target_dir)
                    
                    
                    % take care of magnitude images
                    if strcmp(fmap_type_ls(ifmap_type_ls), 'phasediff') ...
                            || strcmp(fmap_type_ls(ifmap_type_ls), 'phase12')
                        
                        % there might be 2 magnitude images
                        for imag = 1:2
                            
                            [~,name,ext] = spm_fileparts(spm_select('List',path,...
                                sprintf('^.*magnitude%i.*$', imag) ) );
                            
                            in = fullfile(path, [name ext]);
                            
                            % mostly in case there is only one
                            if ~isempty(name)
                                
                                [ext,name,compressed] = iscompressed(ext,name);
                                
                                if imag==1
                                    subjects{s}.fieldmap(ifmap,:).mag1 = ...
                                        fullfile(target_dir, [name ext]);
                                    file_exists = exist(subjects{s}.fieldmap(ifmap,:).mag1,'file');
                                elseif imag==2
                                    subjects{s}.fieldmap(ifmap,:).mag2 = ...
                                        fullfile(target_dir, [name ext]);
                                    file_exists = exist(subjects{s}.fieldmap(ifmap,:).mag2,'file');
                                end
                                
                                unzip_or_copy(compressed, options, file_exists, in, target_dir)
                                
                            end
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
end

return

disp('spmup has finished unpacking data')

%% run preprocessing using options
% -------------------------------------------------------------------------
if isfield(options,'removeNvol') %% need to check BIDS spec (???)
    start_at = options.removeNvol;
else
    start_at = 1;
end

% do the computation if needed
for s=1%:size(subjs_ls,2)
    %     parfor s=1%:size(subjs_ls,2)
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

function [ext,name,compressed] = iscompressed(ext,name)
if strcmp(ext,'.gz')
    compressed = 1;
    [~,name,ext] = spm_fileparts(name);
elseif strcmp(ext,'.nii')
    compressed = 0;
else
    error('Unknown extension type.')
end
end

function unzip_or_copy(compressed, options, file_exists, in, target_dir)
if compressed % if compressed
    if strcmp(options.overwrite_data,'on') || ( strcmp(options.overwrite_data,'off') ...
            && ~file_exists )
        gunzip(in, target_dir);
    end
else % if not compressed
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~file_exists )
        copyfile(in, target_dir);
    end
end
end