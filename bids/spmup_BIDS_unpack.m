function [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options)

% routine to read and unpack BIDS fMRI data
% can also be used to list files to be preprocessed to pass to spmup_BIDS_preprocess
%
% FORMAT spmup_BIDS_unpack
%        spmup_BIDS_unpack(BIDS_dir)
%        spmup_BIDS_unpack(BIDS_dir,options)
%
% INPUT - BIDS_dir is the BIDS directory
%       - options are given using 'key', 'value' pairs in a structure
%         see spmup_getoptions - the minimal required field is 
%         'outdir' where to write the data (default is [BIDS_dir filesep derivatives])
%         additional fields can be used to filter data
%         'task' to choose only some fMRI tasks with value the task name
%         'anat' to choose only some anat files with values as a cell array
%                of pattern e.g. {'rec-distcorr_T1w','rec-distcorr_T2w'} 
%
% OUTPUT - BIDS: the structure returned by spm_BIDS and possibly modified by spmup_BIDS_unpack
%        - subjects: a structure containing the fullpath of the unpacked
%                  anat, fmap and func files for each subject
%
% Example usage BIDS_dir        = 'F:\WakemanHenson_Faces\fmri';
%               options         = spmup_getoptions(BIDS_dir); 
%               options.anat    = {'_T1w'}; % will only unpack the sub-*_T1w.nii images
%               options.task    = {'facerecognition','plt'}; % will only unpack those 2 tasks
%               [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options)
%
%               % you already unpacked but just want the lists
%               BIDS_dir        = 'F:\WakemanHenson_Faces\fmri';
%               options         = spmup_getoptions(BIDS_dir); 
%               options.anat    = {'_T1w'}; % will only lists the sub-*_T1w.nii images
%               options.task    = {'facerecognition'}; % will only lists this task
%               options.overwrite_data = 'off'; % magic !! 
%               [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options)
%
% see also spmup_getoptions
% Cyril Pernet & Remi Gau
% --------------------------
%  Copyright (C) SPMUP Team 

%% inputs
% -------------------------------------------------------------------------
if nargin == 0
    BIDS_dir = uigetdir(pwd,'select BIDS directory');
    if BIDS_dir == 0
        return
    end
end

% update options with user input
% --------------------------------
if ~exist('options','var')
    options = spmup_getoptions(BIDS_dir);
end

if ~isfield(options,'outdir')
    options.outdir = [BIDS_dir filesep 'derivatives'];
end

if ~isfield(options,'task')
    options.task = []; % = do them all
end

% get the data info
% --------------------
spm('defaults', 'FMRI');
disp('getting BIDS info')
BIDS = spm_BIDS(BIDS_dir);

%% unpack data
% -------------------------------------------------------------------------
if ~exist(options.outdir,'dir')
    mkdir(options.outdir);
end

subjs_ls = spm_BIDS(BIDS,'subjects');
nb_sub   = numel(subjs_ls);

% unpack anat and center [0 0 0]
% -------------------------------
for s=nb_sub:-1:1 % for each subject

    sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
    if isempty(sess_ls)
        sess_ls{1} = '';
    end
    BIDSindex = find(arrayfun(@(x) ~isempty(strfind(['sub-' subjs_ls{s}],x.name)), BIDS.subjects));
    
    if ~exist(fullfile(options.outdir, ['sub-' subjs_ls{s}]),'dir')
        mkdir(fullfile(options.outdir, ['sub-' subjs_ls{s}]))
    end

    for session = 1:length(sess_ls) % for each session

        if isempty(sess_ls{session})
            sess_folder = '';
        else
            sess_folder = ['ses-' sess_ls{session}];
        end

        % target directory where the anat images will be copied/unzipped
        target_dir = fullfile(options.outdir, ['sub-' subjs_ls{s}], sess_folder, 'anat'); 
        fprintf(' subject %g - session %i: checking anat \n', s, session)

        % lists all the anat images for that subject and session
        for type = 1:length(options.anat)
            in = char(spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
                'ses', sess_ls{session}, 'modality', 'anat'));
            
            % only pick the anat file matching the pattern
            if ~isempty(options.anat)
                pattern = fullfile(BIDS.subjects(BIDSindex).path, ...
                    ['anat' filesep 'sub-' subjs_ls{s} options.anat{type}]);
                for match = size(in,1):-1:1
                    inindex(match) = strncmp(pattern,in(match,:),length(pattern));
                end
                in = deblank(in(inindex,:)); 
                clear inindex
            end
            
            if isempty(in) % in case there is no anat file for this session
                warning(' No valid anat file found for subject %s - session %s', ...
                    subjs_ls{s}, sess_ls{session})
            else
                for file = 1:size(in,1) % we may want to unpack them all
                    [~,name,ext,~] = spm_fileparts(in(file,:));
                    [ext,name,compressed] = iscompressed(ext,name);
                    subjects{s}.anat{file} = fullfile(target_dir, [name ext]); % keep track of where the file is
                    unzip_or_copy(compressed, options, exist(subjects{s}.anat{file},'file'), in, target_dir)
                end % for each in files
            end 
        end % for each options.anat
    end % for each session
end % for each subject

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
parfor s=1:nb_sub
    sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
    if isempty(sess_ls)
        sess_ls = {''};
    end
    
    if ~exist(fullfile(options.outdir, ['sub-' subjs_ls{s}]),'dir') %#ok<PFBNS>
        mkdir(fullfile(options.outdir, ['sub-' subjs_ls{s}]))
    end

    % counters to list files across different sessions
    bold_run_count = 1;
    fmap_run_count = 1;

    for session = 1:size(sess_ls,2) % for each session

        fprintf(' subject %g - session %i: checking functional data \n',s, session)
        if size(sess_ls,2)==1
            sess_folder = '';
        else
            sess_folder = ['ses-' sess_ls{session}];
        end

        if ~exist(fullfile(options.outdir, ['sub-' subjs_ls{s}], sess_folder, 'func'),'dir')
            mkdir(fullfile(options.outdir, ['sub-' subjs_ls{s}], sess_folder, 'func')); 
        end

        % list all bold files for that subject and session
        run_ls = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
            'ses', sess_ls{session}, 'type', 'bold','task',options.task);
        metadata = spm_BIDS(BIDS, 'metadata', 'sub', subjs_ls{s}, ...
            'ses', sess_ls{session}, 'type', 'bold','task',options.task);

        %% functional
        for frun = 1:size(run_ls,1) % for each run

            target_dir = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                sess_folder, 'func', ['run' num2str(frun)]);            
            in = run_ls{frun,1};
            [~,name,ext] = spm_fileparts(in);
            [ext,name,compressed] = iscompressed(ext,name);
            subjects{s}.func{bold_run_count,1} = fullfile(target_dir, [name ext]);
            if iscell(metadata)
                subjects{s}.func_metadata{bold_run_count,1} = metadata{frun};
            else
                subjects{s}.func_metadata{bold_run_count,1} = metadata;
            end
            unzip_or_copy(compressed, options, exist(subjects{s}.func{bold_run_count,1},'file'), in, target_dir)
            bold_run_count = bold_run_count + 1;
        end

        %% field maps
        if ismember('fmap', spm_BIDS(BIDS,'modalities', 'sub', subjs_ls{s}, ...
                'ses', sess_ls{session}))

            file_exists = 0;
            fprintf(' subject %g - session %i: checking field maps \n',s, session)
            target_dir = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                sess_folder, 'fieldmaps');

            % check which types of fieldmaps we are dealing with
            fmap_type_ls = spm_BIDS(BIDS, 'types', 'sub', subjs_ls{s}, ...
                'ses', sess_ls{session}, 'modality', 'fmap');

            % goes through each type in case we have several fieldmap types in that
            % data set
            for ifmap_type_ls = 1:size(fmap_type_ls,1)

                % lists all the fieldmap files of that type as there might
                % be several runs / acq / rec ...
                fmap_ls = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
                    'ses', sess_ls{session}, 'modality', 'fmap', ...
                    'type', fmap_type_ls{ifmap_type_ls});

                metadata = spm_BIDS(BIDS, 'metadata', 'sub', subjs_ls{s}, ...
                    'ses', sess_ls{session}, 'modality', 'fmap', ...
                    'type', fmap_type_ls{ifmap_type_ls});

                % for all the fieldmaps of that type
                for ifmap = 1:size(fmap_ls,1)

                    in = fmap_ls{ifmap,1};
                    [path,name,ext] = spm_fileparts(in);
                    [ext,name,compressed] = iscompressed(ext,name);
                    subjects{s}.fieldmap(fmap_run_count,1).metadata = ...
                        metadata{ifmap};

                    % different behavior depending on fielpmap type
                    switch fmap_type_ls{ifmap_type_ls}

                        case 'phasediff'
                            subjects{s}.fieldmap(fmap_run_count,1).phasediff = ...
                                fullfile(target_dir, [name ext]);
                            subjects{s}.fieldmap(fmap_run_count,1).type = 'phasediff';
                            file_exists = exist(subjects{s}.fieldmap(fmap_run_count,:).phasediff,'file');

                        case 'phase12'
                            subjects{s}.fieldmap(fmap_run_count,1).type = 'phase12';
                            if contains(name,'phase1')
                                subjects{s}.fieldmap(fmap_run_count,1).phase1 = ...
                                    fullfile(target_dir, [name ext]);
                                file_exists = exist(subjects{s}.fieldmap(fmap_run_count,:).phase1,'file');
                            elseif contains(name,'phase2')
                                subjects{s}.fieldmap(fmap_run_count,1).phase2 = ...
                                    fullfile(target_dir, [name ext]);
                                file_exists = exist(subjects{s}.fieldmap(fmap_run_count,:).phase2,'file');
                            end

                        case 'fieldmap'
                            subjects{s}.fieldmap(fmap_run_count,1).type = 'fieldmap';
                            warning(' Fieldmap type of fielmaps not implemented')

                        case 'epi'
                            subjects{s}.fieldmap(fmap_run_count,1).type = 'epi';
                            warning(' EPI type of fielmaps not implemented')

                        otherwise
                            subjects{s}.fieldmap(fmap_run_count,1).type = 'unknown';
                            warning(' %s is an unsupported type of fieldmap', fmap_type_ls(ifmap_type_ls))
                    end
                    unzip_or_copy(compressed, options, file_exists, in, target_dir)

                    % taking care of magnitude images
                    if strcmp(fmap_type_ls(ifmap_type_ls), 'phasediff') ...
                            || strcmp(fmap_type_ls(ifmap_type_ls), 'phase12')

                        diff_name = name;
                        % there might be 2 magnitude images
                        for imag = 1:2
                            loof_for = strrep(diff_name, fmap_type_ls{ifmap_type_ls}, ...
                                sprintf('magnitude%i', imag) );
                            [~,name,ext] = spm_fileparts(spm_select('List',path,...
                                [loof_for '.*$'] ) );
                            in = fullfile(path, [name ext]);

                            % mostly in case there is only one
                            if ~isempty(name)

                                [ext,name,compressed] = iscompressed(ext,name);
                                if imag==1
                                    subjects{s}.fieldmap(fmap_run_count,1).mag1 = ...
                                        fullfile(target_dir, [name ext]);
                                    file_exists = exist(subjects{s}.fieldmap(fmap_run_count,1).mag1,'file');
                                elseif imag==2
                                    subjects{s}.fieldmap(fmap_run_count,1).mag2 = ...
                                        fullfile(target_dir, [name ext]);
                                    file_exists = exist(subjects{s}.fieldmap(fmap_run_count,1).mag2,'file');
                                end
                                unzip_or_copy(compressed, options, file_exists, in, target_dir)
                            end
                        end
                    end
                    fmap_run_count = fmap_run_count + 1;
                end
            end
        end
    end
end

disp('spmup has finished unpacking data')
try delete(gcp('nocreate')); end %#ok<TRYNC>

end

function [ext,name,compressed] = iscompressed(ext,name)
% returns ext, name file and if data are compressed

if strcmp(ext,'.gz')
    compressed   = 1;
    [~,name,ext] = spm_fileparts(name);
elseif strcmp(ext,'.nii')
    compressed = 0;
elseif isempty(name)
    compressed = [];
    warning('No valid file found')
end
end

function unzip_or_copy(compressed, options, file_exists, filein, target_dir)
% compressed is 1 or 0 indicating if one needs to unzip .gz or copy .nii
% options.overwrite_data is 'on' or 'off' indicating to unzip/copy over existing data
% file_exists is 1 or 0 indicating is the file to unzip or copy is already
%             there! if it exists (1) and options.overwrite_data in 'off' nothing is done
% filein is the actual file to unzip or copy
% target_dir is the location to unzip or copy

if strcmp(options.overwrite_data,'on') || ...
        (strcmp(options.overwrite_data,'off') && ~file_exists)
    
    if ~exist(target_dir,'dir')
        mkdir(target_dir);
    end
    
    if compressed % if compressed
        [~,name,ext]   = spm_fileparts(filein);
        fprintf(' Unpacking data: %s\n', fullfile(target_dir, [name ext]))
        gunzip(filein, target_dir);
    else
        fprintf(' Copying data: %s\n',filein)
        copyfile(filein, target_dir);
    end
end
end
