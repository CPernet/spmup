function [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options)

% routine to read and unpack BIDS fMRI data (and anat)
% can also be used to list files to be preprocessed to pass to spmup_BIDS_preprocess
%
% FORMAT [BIDS,subjects] = spmup_BIDS_unpack
%        [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir)
%        [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options)
%
% INPUT - BIDS_dir is the BIDS directory
%       - options are given using 'key', 'value' pairs in a structure
%         see spmup_getoptions - the minimal required field is 
%         'outdir' where to write the data (default is [BIDS_dir filesep derivatives])
%         additional fields can be used to filter data, see options
%         for instance 'subjects' to choose only some subjects
%
% OUTPUT - BIDS: the structure returned by spm_BIDS and possibly modified by spmup_BIDS_unpack
%        - subjects: a structure containing the fullpath of the unpacked
%                  anat, fmap and func files for each subject + link to events
%                  >>> dimension 1 encodes sessions, dimension 2 scans, runs, etc
%
% Example usage BIDS_dir        = 'F:\WakemanHenson_Faces\fmri';
%               options         = spmup_getoptions(BIDS_dir); 
%               options.anat    = {'_T1w'}; % will only unpack the sub-*_T1w.nii images
%               options.task    = {'facerecognition','plt'}; % will only unpack those 2 tasks
%               [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options)
%               % you may need to run system('chmod 755') on linux servers
%                 depending on matlab configuration
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
%
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

% update BIDS for multiple 'overlapping' fmaps
% these can be filtered durnig the analysis
for s=1:size(BIDS.subjects,2)
    location = fullfile(BIDS.subjects(s).path,'fmap');
    if exist(location,'dir')
        contents = dir(fullfile(location,[BIDS.subjects(s).name '*.nii*']));
        if ~isempty(contents) && isempty(BIDS.subjects(s).fmap)
             fmap_index = size(BIDS.subjects(s).fmap,2)+1;
             for file = 1:size(contents,1)
                 name = fullfile(contents(file).folder,contents(file).name);
                 if ~contains(name,'magnitude')
                     json_content = jsondecode(fileread([name(1:strfind(name,'.nii')) 'json']));
                     
                     if contains(name,'acq-')
                         tmp         = name(strfind(name,'acq-')+4:end);
                         acquisition = tmp(1:min(strfind(tmp,'_'))-1);
                     else
                         acquisition = json_content.ScanningSequence;
                     end
                     
                     if contains(name,'dir-')
                         tmp         = name(strfind(name,'dir-')+4:end);
                         direction   = tmp(1:min(strfind(tmp,'_'))-1);
                     else
                         direction   = json_content.PhaseEncodingDirection;
                     end
                     
                     if contains(name,'run-')
                         tmp         = name(strfind(name,'run-')+4:end);
                         run_number  = tmp(1:min(strfind(tmp,'_'))-1);
                     else
                         run_number  = '';
                     end
                     
                     if contains(name,'_epi')
                         type = 'epi';
                     elseif contains(name,'_phasediff')
                         type = 'phasediff';
                         ext  = name(strfind(name,'_phasediff')+length('_phasediff'):end);
                         if exist([name(1:strfind(name,'_phasediff')) 'magnitude' ext],'file')
                             magnitude = [name(1:strfind(name,'_phasediff')) 'magnitude' ext];
                         elseif exist([name(1:strfind(name,'_phasediff')) 'magnitude1' ext],'file')
                             if exist([name(1:strfind(name,'_phasediff')) 'magnitude2' ext],'file')
                                 magnitude{1} = [name(1:strfind(name,'_phasediff')) 'magnitude1' ext];
                                 magnitude{2} = [name(1:strfind(name,'_phasediff')) 'magnitude2' ext];
                             else
                                 magnitude = [name(1:strfind(name,'_phasediff')) 'magnitude1' ext];
                             end
                         end
                     elseif any(contains(name,{'_phase1','_phase2'}))
                         type = 'phase';
                         tmp  = name; clear name
                         ext  = tmp(strfind(tmp,'_phase')+length('_phase')+1:end);
                         if contains(tmp,'_phase1')
                             if exist([tmp(1:strfind(tmp,'_phase1')) 'phase2' ext],'file')
                                 name{1} = tmp;
                                 name{2} = [tmp(1:strfind(tmp,'_phase1')) 'phase2' ext];
                             else
                                 error('phase1 fmap found but not phase2')
                             end
                         else
                             if exist([tmp(1:strfind(tmp,'_phase2')) 'phase1' ext],'file')
                                 name{1} = [tmp(1:strfind(tmp,'_phase2')) 'phase1' ext];
                                 name{2} = tmp;
                             else
                                 error('phase2 fmap found but not phase1')
                             end
                         end
                         if exist([tmp(1:strfind(tmp,'_phase')) 'magnitude1' ext],'file')
                             if exist([tmp(1:strfind(tmp,'_phase')) 'magnitude2' ext],'file')
                                 magnitude{1} = [tmp(1:strfind(tmp,'_phase')) 'magnitude1' ext];
                                 magnitude{2} = [tmp(1:strfind(tmp,'_phase')) 'magnitude2' ext];
                             else
                                 magnitude = [tmp(1:strfind(tmp,'_phase')) 'magnitude1' ext];
                             end
                         end
                     end
                     
                     BIDS.subjects(s).fmap(fmap_index).type          = type;
                     BIDS.subjects(s).fmap(fmap_index).filename      =  name;
                     BIDS.subjects(s).fmap(fmap_index).ses           =  BIDS.subjects(s).session;
                     BIDS.subjects(s).fmap(fmap_index).acq           =  acquisition;
                     BIDS.subjects(s).fmap(fmap_index).dir           =  direction;
                     BIDS.subjects(s).fmap(fmap_index).run           =  run_number;
                     BIDS.subjects(s).fmap(fmap_index).meta          =  json_content;
                     if ~strcmp(type,'epi')
                         BIDS.subjects(s).fmap(fmap_index).magnitude = magnitude;
                     end
                     fmap_index = fmap_index+1;
                 end
             end
                 
         elseif ~isempty(contents) && ~isempty(BIDS.subjects(s).fmap)
            fmap_index = size(BIDS.subjects(s).fmap,2)+1;
            for file = 1:size(contents,1)
                existing_names = arrayfun(@(x) x.filename, BIDS.subjects(s).fmap, 'UniformOutput',false);
                islisted = any(cell2mat(cellfun(@(x) strcmp(x, contents(file).name), existing_names, 'UniformOutput', false)));
                if ~islisted
                    name = fullfile(contents(file).folder,contents(file).name);
                    if ~contains(name,'magnitude')
                        json_content = jsondecode(fileread([name(1:strfind(name,'.nii')) 'json']));
                        
                        if contains(name,'acq-')
                            tmp         = name(strfind(name,'acq-')+4:end);
                            acquisition = tmp(1:min(strfind(tmp,'_'))-1);
                        else
                            acquisition = json_content.ScanningSequence;
                        end
                        
                        if contains(name,'dir-')
                            tmp         = name(strfind(name,'dir-')+4:end);
                            direction   = tmp(1:min(strfind(tmp,'_'))-1);
                        else
                            direction   = json_content.PhaseEncodingDirection;
                        end
                        
                        if contains(name,'run-')
                            tmp         = name(strfind(name,'run-')+4:end);
                            run_number  = tmp(1:min(strfind(tmp,'_'))-1);
                        else
                            run_number  = '';
                        end
                        
                        if contains(name,'_epi')
                            type = 'epi';
                        elseif contains(name,'_phasediff')
                            type = 'phasediff';
                            ext  = name(strfind(name,'_phasediff')+length('_phasediff'):end);
                            if exist([name(1:strfind(name,'_phasediff')) 'magnitude' ext],'file')
                                magnitude = [name(1:strfind(name,'_phasediff')) 'magnitude' ext];
                            elseif exist([name(1:strfind(name,'_phasediff')) 'magnitude1' ext],'file') 
                                if exist([name(1:strfind(name,'_phasediff')) 'magnitude2' ext],'file') 
                                    magnitude{1} = [name(1:strfind(name,'_phasediff')) 'magnitude1' ext];
                                    magnitude{2} = [name(1:strfind(name,'_phasediff')) 'magnitude2' ext];
                                else
                                    magnitude = [name(1:strfind(name,'_phasediff')) 'magnitude1' ext];
                                end
                            end
                        elseif any(contains(name,{'_phase1','_phase2'}))
                            type = 'phase';
                            tmp  = name; clear name
                            ext  = tmp(strfind(tmp,'_phase')+length('_phase')+1:end);
                            if contains(tmp,'_phase1')
                                if exist([tmp(1:strfind(tmp,'_phase1')) 'phase2' ext],'file')
                                    name{1} = tmp;
                                    name{2} = [tmp(1:strfind(tmp,'_phase1')) 'phase2' ext];
                                else
                                   error('phase1 fmap found but not phase2') 
                                end
                            else
                                if exist([tmp(1:strfind(tmp,'_phase2')) 'phase1' ext],'file')
                                    name{1} = [tmp(1:strfind(tmp,'_phase2')) 'phase1' ext];
                                    name{2} = tmp;
                                else
                                   error('phase2 fmap found but not phase1') 
                                end
                            end
                            if exist([tmp(1:strfind(tmp,'_phase')) 'magnitude1' ext],'file')
                                if exist([tmp(1:strfind(tmp,'_phase')) 'magnitude2' ext],'file')
                                    magnitude{1} = [tmp(1:strfind(tmp,'_phase')) 'magnitude1' ext];
                                    magnitude{2} = [tmp(1:strfind(tmp,'_phase')) 'magnitude2' ext];
                                else
                                    magnitude = [tmp(1:strfind(tmp,'_phase')) 'magnitude1' ext];
                                end
                            end
                        end
                        
                        BIDS.subjects(s).fmap(fmap_index) = struct(...
                            'type',type, ...
                            'filename', name, ...
                            'ses', BIDS.subjects(s).session, ...
                            'acq', acquisition, ...
                            'dir', direction, ...
                            'run', run_number, ...
                            'meta', json_content);
                        if ~strcmp(type,'epi')
                            BIDS.subjects(s).fmap(fmap_index).magnitude = magnitude;
                        end
                        fmap_index = fmap_index+1;
                    end
                end
            end
        end
    end
end

% remove subjects from BIDS
if ~isempty(options.subjects)
    subindex = 1;
    subnames = arrayfun(@(x) x.name,BIDS.subjects,'UniformOutput',false);
    for sub = 1:length(options.subjects)
        if strncmp(options.subjects{sub},'sub-',4)
            tmp = find(cellfun(@(x) strcmpi(x,options.subjects{sub}),subnames));
        else
            tmp = find(cellfun(@(x) strcmpi(x,['sub-' options.subjects{sub}]),subnames));
        end
        
        if ~isempty(tmp)
            keep{subindex} = tmp; %#ok<AGROW>
            subindex = subindex+1;
        end
    end
    
    if exist('keep','var')
        subindex = 1;
        allsub = zeros(sum(cellfun(@(x) length(x),keep)),1);
        for k=1:length(keep)
            allsub(subindex:subindex+length(keep{k})-1) = keep{k};
            subindex = subindex+length(keep{k});
        end
        tmp           = BIDS;
        BIDS          = rmfield(BIDS,'subjects');
        BIDS.subjects = tmp.subjects(allsub);
        clear tmp keep
    
        % also update BIDS.participants
        subindex = 1;
        subnames = arrayfun(@(x) x.participant_id,BIDS.participants,'UniformOutput',false);
        for sub = 1:length(options.subjects)
            if strncmp(options.subjects{sub},'sub-',4)
                tmp = find(cellfun(@(x) strcmpi(x,options.subjects{sub}),subnames{1}));
            else
                tmp = find(cellfun(@(x) strcmpi(x,['sub-' options.subjects{sub}]),subnames{1}));
            end
            
            if ~isempty(tmp)
                keep(subindex) = tmp; %#ok<AGROW>
                subindex = subindex+1;
            end
        end
        
        if exist('keep','var')
            fn = fieldnames(BIDS.participants);
            for f=1:size(fn,1)
                if ~strcmpi(fn{f},'meta')
                    BIDS.participants.(fn{f}) = BIDS.participants.(fn{f})(keep);
                end
            end
        end
        
    else
        error('options.subjects %s was used to filter data, but no match was found')
    end
end


%% unpack data
% -------------------------------------------------------------------------
if ~exist(options.outdir,'dir')
    mkdir(options.outdir);
end

subjs_ls  = spm_BIDS(BIDS,'subjects');
nb_sub    = numel(subjs_ls);

% unpack anat 
% ------------
for s = nb_sub:-1:1 % for each subject

    sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
    if ~isempty(options.ses)
        sess_ls = sess_ls(strcmpi(sess_ls,options.ses));
        if isempty(sess_ls)
            warning('%s session not found, unpacking all sessions',options.ses)
            sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
        end
    end
    
    if isempty(sess_ls)
        sess_ls{1} = '';
    end

    % BIDSindex = find(arrayfun(@(x) ~isempty(strfind(['sub-' subjs_ls{s}],x.name)), BIDS.subjects));
    
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
        fprintf(' subject %g - session %s: checking anat \n', s, sess_ls{session})

        % lists all the anat images for that subject and session
        in = char(spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
            'ses', sess_ls{session}, 'modality', 'anat'));
        
        if isempty(in) % in case there is no anat file for this session
            warning(' No valid anat file found for subject %s - session %s', ...
                subjs_ls{s}, sess_ls{session})
        else
            % only pick the anat file matching the pattern
            if ~isempty(options.anat)
                matchindex = size(in,1);
                for match = size(in,1):-1:1
                    for o=length(options.anat):-1:1
                        test(o) = contains(in(match,:),options.anat{o});
                    end
                    inindex(matchindex) = any(test);
                    matchindex = matchindex-1;
                end
                in = deblank(in(inindex,:));
                clear inindex
            end
            
            % check integrity of list in
            for i=1:size(in,1)
                names{i} = in(i,:);
            end
            in = unique(names)'; 
            clear names
            
            % unpack
            anat_count = 1;
            for file = 1:size(in,1)
                [~,name,ext,~]        = spm_fileparts(in{file});
                [ext,name,compressed] = iscompressed(deblank(ext),name);
                if ~isempty(compressed)
                    subjects{s}.anat{session,anat_count,1} = fullfile(target_dir, [name deblank(ext)]); % keep track of where the file is
                    unzip_or_copy(compressed, options, exist(subjects{s}.anat{file},'file'), in{file}, target_dir)
                    anat_count = anat_count + 1;
                end
            end % for each in files
            clear in
        end
    end % for each session
end % for each subject

% if more than 2 functional run, try to unpack data using multiple cores
% ----------------------------------------------------------------------
if size(spm_BIDS(BIDS, 'runs', 'sub', subjs_ls{s}),2) > 2 %#ok<*UNRCH>
    spmup_setparallel(options.Ncores)
end

% unpack functional and field maps (still need to work out sessions here)
% -----------------------------------------------------------------------
parfor s=1:nb_sub
    sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
    if ~isempty(options.ses)
        sess_ls = sess_ls(strcmpi(sess_ls,options.ses));
        if isempty(sess_ls)
            warning('%s session not found, unpacking all sessions',options.ses)
            sess_ls = spm_BIDS(BIDS,'sessions', 'sub', subjs_ls{s});
        end
    end
    
    if isempty(sess_ls)
        sess_ls = {''};
    end
    
    if ~exist(fullfile(options.outdir, ['sub-' subjs_ls{s}]),'dir') 
        mkdir(fullfile(options.outdir, ['sub-' subjs_ls{s}]))
    end

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
        if isempty(options.task)
            run_ls = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
                'ses', sess_ls{session}, 'type', 'bold');
            metadata = spm_BIDS(BIDS, 'metadata', 'sub', subjs_ls{s}, ...
                'ses', sess_ls{session}, 'type', 'bold');
        else
            run_ls = spm_BIDS(BIDS, 'data', 'sub', subjs_ls{s}, ...
                'ses', sess_ls{session}, 'type', 'bold', 'task', options.task);
            metadata = spm_BIDS(BIDS, 'metadata', 'sub', subjs_ls{s}, ...
                'ses', sess_ls{session}, 'type', 'bold', 'task', options.task);
        end        

        % counters to list files across different sessions
        bold_run_count = 1;
        fmap_run_count = 1;

        %% functional
        for frun = 1:size(run_ls,1) % for each run
            in                                        = run_ls{frun,1};
            [rootpath,name,ext]                       = spm_fileparts(in);
            [ext,name,compressed]                     = iscompressed(ext,name);
            if contains(name,'run-')
                runname    = name(strfind(name,'run-'):end);
                runname    = runname(1:min(strfind(runname,'_'))-1);
                target_dir = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                    sess_folder, 'func', runname);               
            else
                target_dir = fullfile(options.outdir, ['sub-' subjs_ls{s}], ...
                    sess_folder, 'func', ['run' num2str(frun)]);
            end
            subjects{s}.func{session,bold_run_count}  = fullfile(target_dir, [name ext]);
            % values                                    = strfind(name,'_')-1;
            % task                                      = strfind(name,'task-');
            % task                                      = name(task:min(values(values>task)));
            if exist(fullfile(rootpath, [name(1:end-5) '_events.tsv']),'file')
                subjects{s}.event{session,bold_run_count,1} = fullfile(rootpath, [name(1:end-5) '_events.tsv']);
            end
            
            if iscell(metadata)
                subjects{s}.func_metadata{session,bold_run_count} = metadata{frun};
            else
                subjects{s}.func_metadata{session,bold_run_count} = metadata;
            end
            unzip_or_copy(compressed, options, exist(subjects{s}.func{session,bold_run_count},'file'), in, target_dir)
            bold_run_count = bold_run_count + 1;
        end

        %% if exist physio
        
        
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
            for ifmap_type_ls = 1:size(fmap_type_ls,2)

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

                    in                    = fmap_ls{ifmap,1};
                    [path,name,ext]       = spm_fileparts(in);
                    [ext,name,compressed] = iscompressed(ext,name);
                    if iscell(metadata)
                        subjects{s}.fieldmap(session,fmap_run_count).metadata = ...
                        metadata{ifmap};
                    else
                        subjects{s}.fieldmap(session,fmap_run_count).metadata = metadata;
                    end

                    % different behavior depending on fielpmap type
                    switch fmap_type_ls{ifmap_type_ls}

                        case 'phasediff'
                            subjects{s}.fieldmap(session,fmap_run_count).phasediff = ...
                                fullfile(target_dir, [name ext]);
                            subjects{s}.fieldmap(session,fmap_run_count).type = 'phasediff';
                            file_exists = exist(subjects{s}.fieldmap(session,fmap_run_count).phasediff,'file');

                        case 'phase12'
                            subjects{s}.fieldmap(session,fmap_run_count).type = 'phase12';
                            if contains(name,'phase1')
                                subjects{s}.fieldmap(session,fmap_run_count).phase1 = ...
                                    fullfile(target_dir, [name ext]);
                                file_exists = exist(subjects{s}.fieldmap(session,fmap_run_count).phase1,'file');
                            elseif contains(name,'phase2')
                                subjects{s}.fieldmap(fmap_run_count,1).phase2 = ...
                                    fullfile(target_dir, [name ext]);
                                file_exists = exist(subjects{s}.fieldmap(session,fmap_run_count).phase2,'file');
                            end

                        case 'fieldmap'
                            subjects{s}.fieldmap(session,fmap_run_count).type = 'fieldmap';
                            warning(' Fieldmap type of fielmaps not implemented')

                        case 'epi'
                            subjects{s}.fieldmap(session,fmap_run_count).type = 'epi';
                            warning(' EPI type of fielmaps are not SPM supported ...')
                            if contains(name,{'AP','PA','LR','RL'})
                                subjects{s}.fieldmap(session,fmap_run_count).epi = ...
                                    fullfile(target_dir, [name ext]);
                                file_exists = exist(subjects{s}.fieldmap(session,fmap_run_count).epi,'file');
                            end

                        otherwise
                            subjects{s}.fieldmap(session,fmap_run_count).type = 'unknown';
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
                                [loof_for '.nii';loof_for '.*gz']));
                            in = fullfile(path, [name ext]);

                            % mostly in case there is only one
                            if ~isempty(name)

                                [ext,name,compressed] = iscompressed(ext,name);
                                if imag==1
                                    subjects{s}.fieldmap(session,fmap_run_count).mag1 = ...
                                        fullfile(target_dir, [name ext]);
                                    file_exists = exist(subjects{s}.fieldmap(session,fmap_run_count).mag1,'file');
                                elseif imag==2
                                    subjects{s}.fieldmap(session,fmap_run_count).mag2 = ...
                                        fullfile(target_dir, [name ext]);
                                    file_exists = exist(subjects{s}.fieldmap(session,fmap_run_count).mag2,'file');
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

% update the BIDS structure allowing to keep quering using spm_BIDS
sourcedir = BIDS.dir;
for s=1:size(BIDS.subjects,1)
    for r=1:size(BIDS.subjects(s).func,1)
        BIDS.subjects(s).func(r).run = ['run' num2str(r)];
        BIDS.subjects(s).func(r).meta = spm_BIDS(BIDS, 'metadata', ...
            'sub', BIDS.subjects(s).func(r).sub, ...
            'ses', BIDS.subjects(s).func(r).ses, ...
            'task', BIDS.subjects(s).func(r).task, ...
            'type', 'bold');
    end
    BIDS.subjects(s).path = [options.outdir BIDS.subjects(s).path(length(sourcedir)+1:end)];
end
BIDS.dir  = options.outdir;

disp('spmup has finished unpacking data')
try delete(gcp('nocreate')); end %#ok<TRYNC>

end

function [ext,name,compressed] = iscompressed(ext,name)
% returns ext, name file and if data are compressed

compressed = [];
if strcmp(ext,'.gz')
    compressed   = 1;
    [~,name,ext] = spm_fileparts(name);
elseif strcmp(ext,'.nii')
    compressed = 0;
elseif isempty(name)
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
        gunzip(deblank(filein), target_dir);
    else
        fprintf(' Copying data: %s\n',filein)
        if strcmp(options.overwrite_data,'on')
            copyfile(deblank(filein), target_dir, 'f');
        else
            try
                copyfile(deblank(filein), target_dir, 'f');
            end
        end
    end
end
end
