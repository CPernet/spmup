function [BIDS,subjects,options]=spmup_BIDS_unpack(BIDS_dir,varargin)

% routine to read and unpack BIDS fMRI and preprocess data
% can also be used to list files to be preprocessed to pass to spmup_BIDS_preprocess
%
% FORMAT spmup_BIDS_unpack
%        spmup_BIDS_unpack(BIDS_dir)
%        spmup_BIDS_unpack(BIDS_dir,options)
%
% INPUTS
%           - BIDS_dir is the BIDS directory
%
%           - options are given using 'key', 'value' pairs
%               'outdir' where to write the data (default is [BIDS_dir filesep derivatives])
%               'keep_data' is 'off' (default) or 'on' to keep all steps 
%                           off means only the last processed data are available
%               'overwrite_data' 'on' (default) turning it 'off' is useful to restart
%               	             some processing while kepping previous steps
%               'QC' 'on' (default) or 'off' performs quality controls for
%               	 anatomical (T1) scans and EPI time series
%               'removeNvol' number of initial volumes to remove
%               'VDM' cellarray (either a file on drive or a memory mapped variable)
%               'task' [] (default) to specificy which bold task to analyze.
%                   Default is to run them all together. This could be problematic
%                   if the tasks have different acquisiton parameters or dimensions.
%               'rec' [] (default) to specificy which bold reconstruction
%               to analyze. Default is to run them all together.
%               'acq' [] (default) to specificy which bold acquisition to
%               'despike' 'off' (default) or 'on' runs median despiking
%               'drifter' 'off' ('default') or 'on' removes cardiac and respiratory signals using the drifter toolbox
%               'motionexp' 'off' (default) or 'on' compute voltera motion expansion
%               'scrubbing' 'on' (default) or 'off' find outliers in motion derivatives and in globals
%               'norm' 'EPInorm' (default) or 'T1norm' choice of the type of template for normalization
%               'ignore_fieldmaps' 'on' or 'off' (default) to ignore distorsion correction for T1norm
%               'skernel' [8 8 8] by default is the smoothing kernel
%
% OUTPUTS
%           - BIDS: the structure returned by spm_BIDS and possibly modified
%           by spmup_BIDS_unpack
%
%           - subjects: a structure containing the fullpath of the unpacked anat,
%           fmap and func files for each subject
%
%           options: a structure equal to the options made plus the
%           default of all the non specified fields
%
% usage example:
% [BIDS,subjects,options]=spmup_BIDS_unpack(pwd,'removeNvol', 0, 'keep_data', 'off',  ...
%     'overwrite_data', 'on', 'motionexp', 'on', 'scrubbing', 'on', 'norm', 'T1norm', 'skernel', [8 8 8]);
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

%% inputs
% -------------------------------------------------------------------------
if nargin == 0
    BIDS_dir= uigetdir(pwd,'select BIDS directory');
    if BIDS_dir == 0
        return
    end
end

% update options with user input
% --------------------------------
options = get_all_options(BIDS_dir);
for c=1:length(varargin)
    if ischar(varargin{c})
        field = intersect(varargin{c},fieldnames(options));
    else
        field = {};
    end
    
    if ~isempty(field)
        options.(field{1}) = varargin{c+1};
    end
end

% get the data info
% --------------------
spm('defaults', 'FMRI');
disp('getting BIDS info')
BIDS = spm_BIDS(BIDS_dir);
if isfield(options,'ignore_fieldmaps')
    if strcmp(options.ignore_fieldmaps, 'on')
        for s=1:length(BIDS.subjects)
            BIDS.subjects(s).fmap = [];
        end
    end
end


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
nb_sub = numel(subjs_ls);


%% unpack anat and center [0 0 0]
% this assumes that there is no more than one T1W image for each session
% -------------------------------
for s=1:nb_sub % for each subject

    sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
    
    mkdir(fullfile(options.outdir, ['sub-' subjs_ls{s}]))

    for session = 1:size(sess_ls,2) % for each session

        if size(sess_ls,2)==1
            sess_folder = '';
        else
            sess_folder = ['ses-' sess_ls{session}];
        end
        
        mkdir(fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                    sess_folder, 'anat'))

        fprintf('subject %g - session %i: checking anat \n', s, session)

        %  target directory where the anat image will be copied/unzipped
        target_dir = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                    sess_folder, 'anat');

        % lists all the T1w images for that subject and session
        in = char(spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
            'ses', sess_ls{session}, 'type', 'T1w'));

        if isempty(in) % in case there is no anat file for this session
            warning('No valid T1w file found for subject %s - session %s', ...
                subjs_ls{s}, sess_ls{session})
        else
            [~,name,ext,~] = spm_fileparts(in);

            % check if we are dealing with a .nii or .nii.gzz file
            % returns ext='.nii' in any case
            % and filename with no extension
            % compressed is 1 or 0
            [ext,name,compressed] = iscompressed(ext,name);

            % keep track of where the file is
            subjects{s}.anat = fullfile(target_dir, [name ext]); %#ok<*AGROW>
            file_exists = exist(subjects{s}.anat,'file'); % necessary to avoid overwriting

            % unzip or copies the file to its target directory depending on
            % extention / options / file existing
            unzip_or_copy(compressed, options, file_exists, in, target_dir)

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
for s=1:nb_sub
    % parfor s=1:nb_sub

    sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
    
    if isempty(sess_ls)
        sess_ls = {''};
    end

    % bold file counter to list theem across different sessions
    bold_run_count = 1;

    % fmap file counter to list theem across different sessions
    fmap_run_count = 1;

    for session = 1:size(sess_ls,2) % for each session

        fprintf('subject %g - session %i: checking functional data \n',s, session)

        if size(sess_ls,2)==1
            sess_folder = '';
        else
            sess_folder = ['ses-' sess_ls{session}];
        end

        [~,~,~] = mkdir(fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
            sess_folder, 'func'));

        % list all bold files for that subject and session
        run_ls = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
            'ses', sess_ls{session}, 'type', 'bold');

        metadata = spm_BIDS(BIDS, 'metadata', 'sub', subjs_ls{s}, ...
            'ses', sess_ls{session}, 'type', 'bold');

        %% functional
        for frun = 1:size(run_ls,1) % for each run

            target_dir = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                sess_folder, 'func', ['run' num2str(frun)]);
            
            if ~exist(target_dir,'dir')
                mkdir(target_dir);
            end
            
            in = run_ls{frun,1};

            [~,name,ext] = spm_fileparts(in);

            [ext,name,compressed] = iscompressed(ext,name);

            % we keep track of where the files are stored
            subjects{s}.func{bold_run_count,1} = fullfile(target_dir, [name ext]);
            if iscell(metadata)
                subjects{s}.func_metadata{bold_run_count,1} = metadata{frun};
            else
                subjects{s}.func_metadata{bold_run_count,1} = metadata;
            end
            file_exists = exist(subjects{s}.func{bold_run_count,1},'file');

            unzip_or_copy(compressed, options, file_exists, in, target_dir)

            bold_run_count = bold_run_count + 1;
        end


        %% field maps
        if ismember('fmap', spm_BIDS(BIDS,'modalities', 'sub', subjs_ls{s}, ...
                'ses', sess_ls{session}))

            fprintf('subject %g - session %i: checking field maps \n',s, session)

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
                            warning('Fieldmap type of fielmaps not implemented')

                        case 'epi'
                            subjects{s}.fieldmap(fmap_run_count,1).type = 'epi';
                            warning('EPI type of fielmaps not implemented')

                        otherwise
                            subjects{s}.fieldmap(fmap_run_count,1).type = 'unknown';
                            warning('%s is an unsupported type of fieldmap', fmap_type_ls(ifmap_type_ls))
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

try delete(gcp('nocreate'));
catch
end

end


%% option sub-functions
% -------------------------------------------------------------------------
function options = get_all_options(BIDS_dir)
% returns all available options set to run the minimal pipeline

options.outdir           = [BIDS_dir filesep 'derivatives' filesep  'spmup_processed'];
options.keep_data        = 'off' ;
options.overwrite_data   = 'on' ;
options.QC               = 'on';
options.removeNvol       = 0;
options.VDM              = {};
options.task             = {};
options.rec              = {};
options.acq              = {};
options.despike          = 'off';
options.drifter          = 'off';
options.motionexp        = 'off';
options.scrubbing        = 'on';
options.norm             = 'EPInorm';
options.ignore_fieldmaps = 'off' ;
options.skernel          = [8 8 8];
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

function unzip_or_copy(compressed, options, file_exists, in, target_dir)

if isempty(in) % if no valid file was found

else

  if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
          && ~file_exists )
      
      fprintf(' Copying data: %s\n',in)
      if ~exist(target_dir,'dir')
          mkdir(target_dir);
      end
      copyfile(in, target_dir);

      if compressed % if compressed

        [~,name,ext] = spm_fileparts(in);
        file_to_unpack = fullfile(target_dir, [name ext]);
        fprintf(' Unpacking data: %s\n', file_to_unpack)
        gunzip(file_to_unpack, target_dir);

      end

  end

end

end
