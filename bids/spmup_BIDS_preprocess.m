function [subject,options] = spmup_BIDS_preprocess(subject, options)

% routine to preprocess BIDS fMRI data following 'options'
% Preprocessing is intented FOR ONE TASK AND ONE SESSION (but multiple runs ok)
%
% FORMAT spmup_BIDS_preprocess(subject, s)
%        spmup_BIDS_preprocess(subjects s, options)
%
% INPUT - subject: a structure containing the fullpath of the unpacked anat,
%           fmap and func files for a given subject (see spmup_BIDS_unpack)
%           if multiple tasks and sesions, filter the subject structure to
%           data appears the dimension 1 indicating scans, runs, etc (see e.g. run_spmup_bids)
%       - options the structure of SPMUP options (see spmup_getoptions)
%
% OUTPUT subject is the same structure updated
%        - with the preprocessed data information
%        - with QC metrics
%        options, if field maps needs computed the .VDM field is updated
%
% Example usage: BIDS_dir          = 'F:\WakemanHenson_Faces\fmri';
%                options           = spmup_getoptions(BIDS_dir);
%                options.overwrite_data = 'off';
%                [BIDS,subjects]   = spmup_BIDS_unpack(BIDS_dir,options);
%                [subject,options] = spmup_BIDS_preprocess(subjects{7}, options)
%
% Cyril Pernet, Remi Gau & Patrick Fisher
% ---------------------------------------
%  Copyright (C) SPMUP Team

% -------------------
%% check a few things
% -------------------

if strcmpi(options.norm,'T1norm')
    if isempty(options.VDM) && strcmpi(options.fmap,'off')
        warning(' Calhoun et al. (2017) showed T1 normalization is not recommended without FieldMaps \n (%s) add VDM map to options \n', ...
            'https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.23737')
    end
end

% each subject build a job structure around matlabbatch
% this is where each subject is analysed using the various options
meta        = [];  % jsom structure to write
opts.indent = ' '; % for spm_jsonwrite
spm_root    = fileparts(which('spm'));
spm_jobman('initcfg');

% ------------------------------------------------
%% reorient anat file to template, apply to others
% ------------------------------------------------

% check T1w is always in position 1
t1index                           = find(cellfun(@(x) contains(x,'T1w'),subject.anat));
for t1 = 1:size(t1index,1)
    anat{t1}                      = subject.anat{t1index(t1)};
end

othersindex                       = find(cellfun(@(x) ~contains(x,'T1w'),subject.anat));
if ~isempty(othersindex)
    for ot = 1:size(othersindex,1)
        anat{end+1}               = subject.anat{othersindex(ot)};
    end
end
subject.anat                      = rmfield(subject,'anat');
subject.anat                      = anat';
clear anat

% process
RMexist = 0;
[anat_dir,filename] = fileparts(subject.anat{1});
if exist(fullfile(anat_dir,[filename(1:end-4) '_desc-preprocessed_anat.json']),'file')
    anatinfo = spm_jsonread(fullfile(anat_dir,[filename(1:end-4) '_desc-preprocessed_anat.json']));
    if isfield(anatinfo,'reorientation_matrix')
        RMexist = 1;
    end
end

if strcmp(options.overwrite_data,'on') || ...
        (strcmp(options.overwrite_data,'off') && ~RMexist)
    
    Data = [subject.anat; subject.func];
    if isfield(subject,'fieldmap')
        if size(subject.fieldmap,1) == 1
            FMindex = structfun(@(x) ischar(x), subject.fieldmap);
            FN      = fieldnames(subject.fieldmap);
            FN      = FN(FMindex);
            for f= 1:size(FN,1)
                if exist(subject.fieldmap.(FN{f}),'file')
                    Data = [Data ; subject.fieldmap.(FN{f})];
                end
            end
        else % 1 field map per run
            for f= 1:size(subject.fieldmap,1)
                FMindex = structfun(@(x) ischar(x), subject.fieldmap(f));
                FN      = fieldnames(subject.fieldmap(f));
                FN      = FN(FMindex);
                for ff= 1:size(FN,1)
                    if exist(subject.fieldmap(f).(FN{ff}),'file')
                        Data = [Data ; subject.fieldmap(f).(FN{ff})];
                    end
                end
            end
        end
    end
    
    if ~isempty(options.VDM)
        for v = 1:length(options.VDM)
            Data = [Data ; options.VDM{v}];
        end
    end
    
    if ~RMexist
        RM = spmup_auto_reorient(Data,1);
        disp(' Data reoriented');
    else
        warning('options.overwrite_data is %s,\nbut data are already reoriented = skipping this step',options.overwrite_data)
    end
    
    % update info about all those files
    func_index = 0;
    fmap_index = 0;
    for f=1:size(Data,1)
        [filepath,filename] = fileparts(Data{f});
        if f==1 && ~RMexist
            anatinfo.reoriented           = '(0,0,0) set with spmup_auto_reorient';
            anatinfo.reorientation_matrix = RM;
            spm_jsonwrite(fullfile(filepath,[filename(1:end-4) '_desc-preprocessed_anat.json']),anatinfo,opts)
        else
            if contains(filepath,'anat')
                meta.reoriented = '(0,0,0) set with spmup_auto_reorient';
                spm_jsonwrite(fullfile(filepath,[filename '_desc-preprocessed_anat.json']),meta,opts)
            elseif contains(filepath,'func')
                if func_index == 0
                    func_index = 1;
                end
                meta.reoriented = '(0,0,0) set from anat with spmup_auto_reorient';
                spm_jsonwrite([subject.func{func_index}(1:end-9) '_desc-preprocessed_bold.json'],meta,opts)
                func_index = func_index+1;
            elseif contains(filepath,'fmap')
                if fmap_index == 0
                    fmap_index = 1;
                end
                meta.reoriented = '(0,0,0) set from anat with spmup_auto_reorient';
                spm_jsonwrite([subject.fmap{fmap_index}(1:end-9) '_desc-preprocessed-fmap.json'],meta,opts)
                fmap_index = fmap_index+1;
            end
        end
        clear meta
    end
end

% ----------------------------------------
%% Despiking and slice timing for each run
% ----------------------------------------

bold_include = [];
included_idx = 0;
for frun = 1:size(subject.func, 1) % each run
    
    filesin                 = subject.func{frun};
    [filepath,filename,ext] = fileparts(filesin);
    if exist(fullfile(filepath,[filename '.json']),'file')
        meta = spm_jsonread(fullfile(filepath,[filename '.json']));
    elseif exist(fullfile(filepath,[filename '_space-IXI549_desc-preprocessed_bold.json']),'file')
        meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']));
    elseif exist(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']),'file')
        meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']));
    end
    meta.run      = frun;
    if isfield(subject.func_metadata{frun},'TaskName')
        meta.TaskName = subject.func_metadata{frun}.TaskName;
    else
        start         = strfind(filename,'task');
        keyindices    = strfind(filename,'_');
        keyindices    = min(keyindices(keyindices>start));
        task          = filename(start+length('task-'):keyindices-1);
        subject.func_metadata{frun}.TaskName = task;
        meta.TaskName = subject.func_metadata{frun}.TaskName;
    end
    
    % check that this BOLD file is of the right task, acquisition,
    % reconstruction
    if ~isempty(options.acq)
        acq = contains(filename, ['acq-' options.acq]);
    else
        acq = 1;
    end
    
    if ~isempty(options.rec)
        rec = contains(filename, ['rec-' options.rec]);
    else
        rec = 1;
    end
    
    if ~isempty(options.task)
        task = contains(filename, ['task-' options.task]);
    else
        task = 1;
    end
    
    % only preprocess this file if it fits what has been requested
    bold_include(frun) = all([task acq rec]);
    if bold_include(frun)
        
        included_idx = included_idx + 1;
        if ~isfield(subject.func_metadata{frun}, 'SliceTiming')
            error('Slice Timing Information missing')
        else
            SliceTiming = subject.func_metadata{frun}.SliceTiming;
        end
        
        if ~isfield(subject.func_metadata{frun}, 'RepetitionTime')
            error('RepetitionTime Information missing')
        else
            RepetitionTime = subject.func_metadata{frun}.RepetitionTime;
        end
        
        % -----------------------------------------
        % remove first volumes
        % -----------------------------------------
        
        if options.removeNvol ~=0
            if strcmp(options.overwrite_data,'on')...
                    || ((strcmp(options.overwrite_data,'off') && ~isfield(meta,'NumberOfVolumesDiscarded')))
                
                % remove from the 4D the files we don't want and proceed
                fprintf('\nadjusting 4D file size run %g \n',frun)
                fmri_out                = spmup_skip(filesin,options.removeNvol);
                [filepath,filename,ext] = fileparts(fmri_out);
                
                % write a json file containing the details of what volumes were
                % removed (see BIDS derivatives specs)
                meta.NumberOfVolumesDiscarded = options.removeNvol;
                spm_jsonwrite([subject.func{frun}(1:end-9) '_desc-preprocessed_bold.json'],meta,opts)
            end
        end
        
        % -----------------------------------------
        % despiking using adaptive median filter
        % -----------------------------------------
        
        if strcmp(options.despike,'before')
            if strcmp(options.overwrite_data,'on') || ...
                    (strcmp(options.overwrite_data,'off') && ~isfield(meta,'despike'))
                
                flags              = struct('auto_mask','on', 'method','median', 'window', [],'skip',0, 'savelog', 'off');
                [V,~,loginfo]      = spmup_despike(fullfile(filepath,[filename,ext]),[],flags);
                filesin            = V.fname; clear V;
                meta.Despiked      = 'before slice timing using spmup_despike';
                meta.despike.param = 'median filter based on the autocorrelation function';
                meta.despike.RMS   = loginfo.RMS;
                spm_jsonwrite([subject.func{frun}(1:end-9) '_desc-preprocessed_bold.json'],meta,opts)
            else
                if exist(fullfile(filepath, [filename(1:end-5) '_rec-despiked_bold' ext]),'file')
                    filesin = fullfile(filepath, [filename(1:end-5) '_rec-despiked_bold' ext]);
                end
            end
        end
        
        % -------------
        % slice timing
        % -------------
        % sliceorder - vector containig the acquisition time for each slice in milliseconds
        % refslice   - time in milliseconds for the reference slice
        % timing     - [0 TR]
        
        if exist('SliceTiming', 'var')
            sliceorder  = SliceTiming; % time
            refslice    = sliceorder(round(length(SliceTiming)/2));
            timing      = [0 RepetitionTime];
            
            [filepath,filename,ext]  = fileparts(filesin);
            st_files{included_idx,1} = fullfile(filepath, ['st_' filename ext]); %#ok<*AGROW>
            
            if strcmp(options.overwrite_data, 'on') || ...
                    (strcmp(options.overwrite_data, 'off') && ~isfield(meta,'slicetiming'))
                fprintf(' starting slice timing correction run %g \n',frun)
                spm_slice_timing(filesin, sliceorder, refslice, timing, 'st_');
                meta.slicetiming.reference = refslice;
                spm_jsonwrite([subject.func{frun}(1:end-9) '_desc-preprocessed_bold.json'],meta,opts)
            end
        else
            st_files{included_idx,1} = fullfile(filepath, [filename ext]);
        end
    end
    clear meta
end % end processing per run

% ------------------------
%% Field Map - compute VDM
% ------------------------

if isfield(subject, 'fieldmap') && strcmpi(options.fmap,'on')
    for ifmap = 1:numel(subject.fieldmap)
        
        % index which fieldmap goes to which bold run and vice versa
        if length(subject.fieldmap) == 1 % 1 field map for all
            if isfield(subject.fieldmap(ifmap).metadata,'IntendedFor')
                file_to_parse = subject.fieldmap(ifmap).metadata.IntendedFor(ifmap);
                if iscell(file_to_parse)
                    file_to_parse = file_to_parse{1};
                end
                [~, filename]                                      = spm_fileparts(file_to_parse);
                which_func_file                                    = contains(subject.func, filename);
                subject.fieldmap(ifmap).metadata.which_func(ifmap) = find(which_func_file); % to know which func files needs this fmap
                subject.which_fmap(find(which_func_file))          = ifmap; % to know which fmap is needed by which func files
            else
                which_func_file                                    = ifmap;
                subject.fieldmap(ifmap).metadata.which_func        = ifmap;
                subject.which_fmap(ifmap)                          = ifmap;
            end
        else
            if isfield(subject.fieldmap(ifmap).metadata,'IntendedFor')
                file_to_parse = subject.fieldmap(ifmap).metadata.IntendedFor(ifmap);
                if iscell(file_to_parse)
                    file_to_parse = file_to_parse{1};
                end
                [~, filename]                                      = spm_fileparts(file_to_parse);
                which_func_file                                = contains(subject.func, filename);
                subject.fieldmap(ifmap).metadata.which_func    = find(which_func_file);
                subject.which_fmap(find(which_func_file))      = ifmap;
            else
                which_func_file                                = ifmap;
                subject.fieldmap(ifmap).metadata.which_func    = ifmap;
                subject.which_fmap(ifmap)                      = ifmap;
            end
        end
        
        % only continue if this fmap is needed for any of the bold run for
        % this task / acq / rec
        if any(bold_include(which_func_file))
            
            % computes average of that run
            % finding bold reference run to which the fielmap will be coregistered
            [filepath,filename,ext]  = fileparts(st_files{find(bold_include, 1)});
            avg = fullfile(filepath,[filename(1:end-5) '-mean_bold' ext]);
            if strcmp(options.overwrite_data,'on') || ...
                    (strcmp(options.overwrite_data,'off') && ~exist(avg,'file'))
                spmup_basics(st_files{find(bold_include, 1)},'mean'); % use the mean slice timed EPI for QC
            end
            
            % VDM creation: output images
            % check which types of fieldmaps we are dealing with
            switch subject.fieldmap(ifmap).type
                case 'phasediff'
                    % two magnitudes (use only 1) and 1 phase difference image
                    [filepath,name,ext]  = fileparts(subject.fieldmap(ifmap).phasediff);
                    options.VDM{ifmap,1} = [filepath filesep 'vdm5_sc' name ext];
                case 'phase12'
                    % two magnitude images and 2 phase images
                    [filepath,name,ext]  = fileparts(subject.fieldmap(ifmap).phase1);
                    options.VDM{ifmap,1} = [filepath filesep 'vdm5_sc' name ext];
                case 'fieldmap'
                    warning('Fieldmap type of fielmaps not implemented')
                case 'epi'
                    warning('EPI type of fielmaps not implemented')
                otherwise
                    warning('%s is an unsupported type of fieldmap', subject.fieldmap(ifmap).type)
            end
            
            
            if strcmp(subject.fieldmap(ifmap).type, 'phasediff') || ...
                    strcmp(subject.fieldmap(ifmap).type, 'phase12')
                
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                        && ~exist(options.VDM{ifmap},'file'))
                    
                    % quick fix: 1 phase difference image no magnitude
                    if strcmp(subject.fieldmap(ifmap).type, 'phasediff') && ~isfield(subject.fieldmap(ifmap),'mag1') || ...
                            strcmp(subject.fieldmap(ifmap).type, 'phasediff') && isempty(subject.fieldmap(ifmap).mag1)
                        clear matlabbatch
                        t1index                                                = find(cellfun(@(x) contains(x,'T1w'),subject.anat));
                        if length(t1index>1)
                            t1index = t1index(ifmap); % this should match session order as spm_BIDS
                        end
                        matlabbatch{1}.spm.spatial.coreg.write.ref             = {subject.fieldmap(ifmap).phasediff};
                        matlabbatch{1}.spm.spatial.coreg.write.source          = {subject.anat{t1index}}; %#ok<FNDSB>
                        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
                        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap   = [0 0 0];
                        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask   = 0;
                        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
                        out                                                    = spm_jobman('run',matlabbatch);
                        destination                                            = fileparts(subject.fieldmap(ifmap).phasediff);
                        [source,filename,ext]                                  = spm_fileparts(out{1}.rfiles{1});
                        subject.fieldmap(ifmap).mag1                           = fullfile(destination,[filename ext]);
                        movefile(fullfile(source,[filename ext]),subject.fieldmap(ifmap).mag1)
                        clear matlabbatch;
                    end
                    
                    % coregister fmap to bold reference
                    matlabbatch{1}.spm.spatial.coreg.estimate.ref                 = {avg};
                    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep      = [16 8 4 2 1];
                    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol      = ...
                        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
                    
                    if strcmp(subject.fieldmap(ifmap).type, 'phasediff')
                        % two magnitudes (use only 1) and 1 phase difference image
                        % possibly the magnitude image is missing !! let's use the t1
                        matlabbatch{end}.spm.spatial.coreg.estimate.source = ...
                            {subject.fieldmap(ifmap).mag1};
                        matlabbatch{end}.spm.spatial.coreg.estimate.other = ...
                            {subject.fieldmap(ifmap).phasediff};
                        
                    elseif strcmp(subject.fieldmap(ifmap).type, 'phase12')
                        % two magnitide images and 2 phase images
                        matlabbatch{end}.spm.spatial.coreg.estimate.source = ...
                            {subject.fieldmap(ifmap).mag1};
                        matlabbatch{end}.spm.spatial.coreg.estimate.other = {...
                            subject.fieldmap(ifmap).mag2;...
                            subject.fieldmap(ifmap).phase1;...
                            subject.fieldmap(ifmap).phase2};
                    end
                    
                    fprintf(' coregistering fmap %g to its reference run\n',ifmap)
                    spm_jobman('run',matlabbatch); clear matlabbatch;
                    
                    % VDM creation: input images
                    
                    % quickly check if we have a magnitude image! if not
                    % reslice T1 and put it in there
                    
                    
                    if strcmp(subject.fieldmap(ifmap).type, 'phasediff')
                        % two magnitudes (use only 1) and 1 phase difference image
                        matlabbatch = spmup_get_FM_workflow('phasediff');
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = ...
                            {subject.fieldmap(ifmap).phasediff};
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = ...
                            {subject.fieldmap(ifmap).mag1};
                    elseif strcmp(subject.fieldmap(ifmap).type, 'phase12')
                        % two magnitide images and 2 phase images
                        matlabbatch = spmup_get_FM_workflow('phase&mag');
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase = ...
                            {subject.fieldmap(ifmap).phase1};
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag   = ...
                            {subject.fieldmap(ifmap).mag1};
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase  = ...
                            {subject.fieldmap(ifmap).phase2};
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag    = ...
                            {subject.fieldmap(ifmap).mag2};
                    end
                    
                    % fieldmap metadata
                    fmap_metadata = subject.fieldmap(ifmap).metadata;
                    echotimes =  1000.*[fmap_metadata.EchoTime1 ...
                        fmap_metadata.EchoTime2]; % change to ms
                    matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = echotimes;
                    
                    % func run metadata: if the fmap is applied to several
                    % runs we take the metadata of the first run it must
                    % applied to
                    func_metadata = subject.func_metadata{which_func_file};
                    if numel(func_metadata)>1
                        func_metadata = func_metadata{1};
                    end
                    
                    if isfield(func_metadata,'TotalReadoutTime')
                        TotalReadoutTime = func_metadata.TotalReadoutTime;
                    elseif isfield(func_metadata,'RepetitionTime')
                        TotalReadoutTime = func_metadata.RepetitionTime;
                    elseif isfield(func_metadata,'EffectiveEchoSpacing')
                        TotalReadoutTime = (func_metadata.NumberOfEchos-1)...
                            *func_metadata.EffectiveEchoSpacing;
                    end
                    
                    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = TotalReadoutTime;
                    
                    if isfield(fmap_metadata,'PulseSequenceType')
                        if sum(strfind(fmap_metadata.PulseSequenceType,'EPI')) ~= 0
                            matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
                        else
                            matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
                        end
                    else
                        disp('Field maps is using default sequence! assuming non-EPI acquisition')
                    end
                    
                    if isfield(fmap_metadata,'PhaseEncodingDirection')
                        if strcmp(fmap_metadata.PhaseEncodingDirection,'j') ...
                                || strcmp(fmap_metadata.PhaseEncodingDirection,'y')
                            matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = 1;
                        elseif strcmp(fmap_metadata.PhaseEncodingDirection,'-j') ...
                                || strcmp(fmap_metadata.PhaseEncodingDirection,'-y')
                            matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = -1;
                        end
                    elseif isfield(func_metadata,'PhaseEncodingDirection')
                        if strcmp(func_metadata.PhaseEncodingDirection,'j') ...
                                || strcmp(func_metadata.PhaseEncodingDirection,'y')
                            matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = 1;
                        elseif strcmp(func_metadata.PhaseEncodingDirection,'-j') ...
                                || strcmp(func_metadata.PhaseEncodingDirection,'-y')
                            matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = -1;
                        end
                    else
                        error('No phase encoding direction found')
                    end
                    
                    matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {avg}; % use the mean despiked / slice timed image
                    
                    fprintf(' computing voxel displacement map %g\n',ifmap)
                    out = spm_jobman('run',matlabbatch); clear matlabbatch;
                end
            end
        end
    end
end

% ------------------------
%% Realignment across runs
% ------------------------

if isempty(options.VDM)
    
    % file out of realign will be
    [filepath,filename,ext] = fileparts(st_files{1});
    mean_realigned_file     = [filepath filesep 'mean' filename ext];
    realign_done            = zeros(size(st_files,1),1);
    for frun = size(st_files,1):-1:1
        realigned_files{frun,1}    = st_files{frun}; % because we don't reslice, simple encode the linear transform in the header
        [filepath,filename]        = fileparts(st_files{frun});
        subject.motionfile{frun,1} = [filepath filesep 'rp_' filename '.txt'];
        [filepath,filename]        = fileparts(subject.func{frun});
        if exist(fullfile(filepath,[filename '.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename '.json']));
        elseif exist(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']));
        elseif exist(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']));
        end
        realign_done(frun) = isfield(meta,'realign');
    end
    
    % actually do it
    if strcmp(options.overwrite_data,'on') || ...
            (strcmp(options.overwrite_data,'off') && sum(realign_done)==0)
        
        fprintf('Starting realignment \n')
        for frun = size(st_files,1):-1:1
            matlabbatch{1}.spm.spatial.realign.estwrite.data{frun} = {st_files{frun}}; %#ok<CCAT1>
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
        out = spm_jobman('run',matlabbatch); clear matlabbatch;
    end
    
else % --------------------------------------------------------------------
    
    % which fieldmap for each run of this task / acq /rec
    if iscell(subject.which_fmap)
        which_fmap          = subject.which_fmap(find(bold_include));
    else
        which_fmap          = subject.which_fmap;
    end
    
    % file out of realign will be
    [filepath,filename,ext] = fileparts(st_files{1});
    mean_realigned_file     = [filepath filesep 'meanur' filename ext];
    realign_done            = zeros(size(st_files,1),1);
    for frun = size(st_files,1):-1:1
        [filepath,filename,ext]    = fileparts(st_files{frun});
        realigned_files{frun,1}    = [filepath filesep 'ur' filename ext]; % because we have the reslice here (not linear)
        subject.motionfile{frun,1} = [filepath filesep 'rp_' filename '.txt'];
        [filepath,filename]        = fileparts(subject.func{frun});
        if exist(fullfile(filepath,[filename '.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename '.json']));
        elseif exist(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']));
        elseif exist(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']));
        end
        realign_done(frun) = isfield(meta,'realign');
    end
    
    % actually do it
    if strcmp(options.overwrite_data,'on') || ...
            (strcmp(options.overwrite_data,'off') && sum(realign_done)==0)
        
        fprintf(' starting realignment and unwarping \n')
        for frun = size(st_files,1):-1:1
            matlabbatch{1}.spm.spatial.realignunwarp.data(frun).scans = {st_files{frun}}; %#ok<CCAT1>
            if ~isempty(options.VDM)
                if length(which_fmap)==1
                    matlabbatch{end}.spm.spatial.realignunwarp.data(frun).pmscan = {options.VDM{which_fmap}}; %#ok<CCAT1>
                else
                    matlabbatch{end}.spm.spatial.realignunwarp.data(frun).pmscan = {options.VDM{which_fmap(frun)}}; %#ok<CCAT1>
                end
            else
                matlabbatch{end}.spm.spatial.realignunwarp.data(frun).pmscan = {''};
            end
        end
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.quality    = 1;
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.sep        = 4;
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.fwhm       = 5;
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.rtm        = 1;
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.einterp    = 4;
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.ewrap      = [0 0 0];
        matlabbatch{end}.spm.spatial.realignunwarp.eoptions.weight     = '';
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.basfcn   = [12 12];
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.lambda   = 100000;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.jm       = 0;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.fot      = [4 5];
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.sot      = [];
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.uwfwhm   = 4;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.rem      = 1;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.noi      = 5;
        matlabbatch{end}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{end}.spm.spatial.realignunwarp.uwroptions.uwwhich  = [2 1];
        matlabbatch{end}.spm.spatial.realignunwarp.uwroptions.rinterp  = 4;
        matlabbatch{end}.spm.spatial.realignunwarp.uwroptions.wrap     = [0 0 0];
        matlabbatch{end}.spm.spatial.realignunwarp.uwroptions.mask     = 1;
        matlabbatch{end}.spm.spatial.realignunwarp.uwroptions.prefix   = 'ur';
        out = spm_jobman('run',matlabbatch); clear matlabbatch;
    end
end


% -----------------------------------------------------------------------
% despiking after rather than before realign using adaptive median filter
% -----------------------------------------------------------------------

if strcmp(options.despike,'after')
    for frun = size(subject.func,1):-1:1
        [filepath,filename] = fileparts(subject.func{frun});
        if exist(fullfile(filepath,[filename '.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename '.json']));
        elseif exist(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']));
        elseif exist(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']));
        end
        
        if strcmp(options.overwrite_data,'on') || ...
                (strcmp(options.overwrite_data,'off') && ~isfield(meta,'despike'))
            realigned = dir(fullfile(filepath,['urst_' filename(1:end-5) '*bold.nii']));
            if isempty(realigned)
                realigned = dir(fullfile(filepath,['st_' filename(1:end-5) '*bold.nii']));
            end
            if ~isempty(realigned)
                realigned = fullfile(realigned.folder,realigned.name);
                flags              = struct('auto_mask','on', 'method','median', 'window', [],'skip',0, 'savelog', 'off');
                [~,~,loginfo]      = spmup_despike(fullfile(realigned.folder,realigned.name),[],flags);
                meta.Despiked      = 'after realignemnt using spmup_despike';
                meta.despike.param = 'median filter based on the autocorrelation function';
                meta.despike.RMS   = loginfo.RMS;
                spm_jsonwrite([subject.func{frun}(1:end-9) '_desc-preprocessed_bold.json'],meta,opts)
            end
        end
    end
end

% update the metadata and do QC
% -----------------------------
for frun = 1:size(subject.func,1)
    
    [filepath,filename] = fileparts(subject.func{frun});
    if exist(fullfile(filepath,[filename '.json']),'file')
        meta = spm_jsonread(fullfile(filepath,[filename '.json']));
    elseif exist(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']),'file')
        meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']));
    elseif exist(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']),'file')
        meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']));
    end
    
    if strcmp(options.QC,'on')
        if exist('out','var')
            if isfield(out{1}.sess(frun),'uwrfiles') % write unwrap realign images
                realigned = out{1}.sess(frun).uwrfiles{1};
            elseif exist(cell2mat(out{1}.sess(frun).rfiles),'file') % write realign images
                realigned = out{1}.sess(frun).rfiles{1};
            else
                realigned = out{1}.sess(frun).cfiles{1}; % realign no writing (not needed = st files)
            end
        else
            [filepath,filename] = fileparts(subject.func{frun});
            realigned = dir(fullfile(filepath,['urst_' filename(1:end-5) '*bold.nii']));
            if isempty(realigned)
                realigned = dir(fullfile(filepath,['st_' filename(1:end-5) '*bold.nii']));
            end
            if ~isempty(realigned)
                realigned = fullfile(realigned.folder,realigned.name);
            else
                error('no realigned file found run %g',frun)
            end
        end
        
        if ~isempty(realigned) % after renaming file, re-running with overwrite 'off' that file is gone
            if ~isfield(subject,'func_qa')
                subject.func_qa                     = cell(size(subject.func));
            end
            if ~isfield(meta,'QA')
                meta.QA.spatialcorr.volume_outliers = cell(size(subject.func));
                meta.QA.spatialcorr.slice_outliers  = cell(size(subject.func));
            end
            
            if isempty(subject.func_qa{frun})
                [subject.func_qa{frun}.volume_outliers, subject.func_qa{frun}.slice_outliers] = spmup_spatialcorr(realigned);
                if iscell(meta.QA.spatialcorr.volume_outliers)
                    meta.QA.spatialcorr.volume_outliers{frun} = subject.func_qa{frun}.volume_outliers;
                    meta.QA.spatialcorr.slice_outliers{frun}  = subject.func_qa{frun}.slice_outliers;
                else
                    meta.QA.spatialcorr.volume_outliers = subject.func_qa{frun}.volume_outliers;
                    meta.QA.spatialcorr.slice_outliers  = subject.func_qa{frun}.slice_outliers;
                end
            else
                if iscell(meta.QA.spatialcorr.volume_outliers)
                    if isempty(meta.QA.spatialcorr.volume_outliers{frun})
                        [meta.QA.spatialcorr.volume_outliers{frun}, meta.QA.spatialcorr.slice_outliers{frun}] = spmup_spatialcorr(realigned);
                    end
                else
                    if isempty(meta.QA.spatialcorr.volume_outliers)
                        [meta.QA.spatialcorr.volume_outliers, meta.QA.spatialcorr.slice_outliers] = spmup_spatialcorr(realigned);
                    end
                end
                subject.func_qa{frun}.volume_outliers = meta.QA.spatialcorr.volume_outliers{frun};
                subject.func_qa{frun}.slice_outliers  = meta.QA.spatialcorr.slice_outliers{frun};
            end
        end
    end
    
    if isempty(options.VDM)
        meta.realign.type       = 'realign without field map correction';
    else
        meta.realign.type       = 'realign and unwarp (with field map correction)';
    end
    spm_jsonwrite([subject.func{frun}(1:end-9) '_desc-preprocessed_bold.json'],meta,opts)
end

% -------------------------------
%% Coregistration / Segmentation 
% --------------------------------

[filepath,filename,ext] = fileparts(subject.anat{1}); % T1 is now always in starting position
spaceEPI_T1w            = fullfile(filepath, ['space-EPI' filename ext]);
filename                = filename(1:strfind(filename,'_T1w')-1);
if ~exist('anatinfo','var') && exist(fullfile(filepath,[filename '_desc-preprocessed_anat.json']),'file')
    anatinfo = spm_jsonread(fullfile(filepath,[filename '_desc-preprocessed_anat.json']));
end

if strcmp(options.overwrite_data,'on') || ...
        (strcmp(options.overwrite_data,'off') && ~isfield(meta,'segment'))
    
    fprintf(' Coregister and segment \n')
    if exist('matlabbatch','var')
        clear matlabbatch
    end
    
    % coregister and reslice other anatomical data to T1
    % --------------------------------------------------
    if size(subject.anat,1) > 1
        for file = 2:size(subject.anat,1)
            matlabbatch{file-1}.spm.spatial.coreg.write.ref          = {subject.anat{1}};
            matlabbatch{end}.spm.spatial.coreg.write.source          = {subject.anat{file}};
            matlabbatch{end}.spm.spatial.coreg.write.roptions.interp = 4;
            matlabbatch{end}.spm.spatial.coreg.write.roptions.wrap   = [0 0 0];
            matlabbatch{end}.spm.spatial.coreg.write.roptions.mask   = 0;
            matlabbatch{end}.spm.spatial.coreg.write.roptions.prefix = 'space-T1w';
            others = spm_jobman('run',matlabbatch); clear matlabbatch;
        end
    end
    
    % coregister anatomical data to mean EPI
    % ---------------------------------
    matlabbatch{1}.spm.spatial.coreg.estimate.ref                 = {mean_realigned_file};
    matlabbatch{end}.spm.spatial.coreg.estimate.source            = {subject.anat{1}};
    if exist('others','var')
        matlabbatch{end}.spm.spatial.coreg.estimate.other         = others{1}.rfiles;
    else
        matlabbatch{end}.spm.spatial.coreg.estimate.other         = {''};
    end
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep      = [16 8 4 2 1];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol      = ...
        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
    
    % reslice anatomical to mean EPI
    % not used in normalization but good for visualization
    matlabbatch{2}.spm.spatial.coreg.write.ref               = {mean_realigned_file};
    matlabbatch{end}.spm.spatial.coreg.write.source          = {subject.anat{1}};
    matlabbatch{end}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{end}.spm.spatial.coreg.write.roptions.wrap   = [0 0 0];
    matlabbatch{end}.spm.spatial.coreg.write.roptions.mask   = 0;
    matlabbatch{end}.spm.spatial.coreg.write.roptions.prefix = 'space-EPI';
    spm_jobman('run',matlabbatch); clear matlabbatch;
    
    anatinfo.coregistration.ref    = mean_realigned_file;
    if exist('others','var')
        anatinfo.coregistration.source = {subject.anat(1) others{1}.rfiles};
    else
        anatinfo.coregistration.source = {subject.anat(1)};
    end
    spm_jsonwrite(fullfile(filepath,[filename '_desc-preprocessed_anat.json']),anatinfo,opts)
end

% ---------------
%% Normalization
% ---------------

[filepath,filename,ext] = fileparts(subject.anat{1});
NormalizedAnat_file     = [filepath filesep 'wm' filename ext];
class{1}                = [filepath filesep 'c1' filename ext];
class{2}                = [filepath filesep 'c2' filename ext];
class{3}                = [filepath filesep 'c3' filename ext];
Normalized_class{1}     = [filepath filesep 'wc1' filename ext];
Normalized_class{2}     = [filepath filesep 'wc2' filename ext];
Normalized_class{3}     = [filepath filesep 'wc3' filename ext];

for frun = size(subject.func,1):-1:1
    [filepath,filename,ext]  = fileparts(realigned_files{frun});
    Normalized_files{frun,1} = [filepath filesep 'w' filename ext];
    if any(strcmpi(options.norm_res ,{'EPI','EPI-iso'}))
        hdr                      = spm_vol(realigned_files{frun});
        norm_res                 = abs(round(diag(hdr(1).mat)));
        norm_res(end)            = [];
        if ~isfield(meta,'normalise')
            if strcmpi(options.norm_res ,'EPI-iso') && length(unique(norm_res)) ~= 1 % i.e. it is not isotropic
                if length(unique(norm_res)) == 3 % all different, oh boy - round to smallest
                    all_res(frun,:) = repmat(min(unique(norm_res)),[3 1]);
                else
                    iso = cell2mat(arrayfun(@(x) norm_res == x, unique(norm_res), 'UniformOutput', false)');
                    [~,position] = max(sum(iso)); % pick resolution that is square already
                    all_res(frun,:) = repmat(unique(norm_res(iso(:,position))),[3 1]);
                end
            else % EPI
                all_res(frun,:) = norm_res;
            end
        else
            all_res(frun,:) = meta.normalise.resolution;
        end
    else % actual value requested
        if ischar(options.norm_res)
            options.norm_res = str2double(options.norm_res);
        end
        if length(options.norm_res) == 1
            options.norm_res = repmat(options.norm_res,[1 3]);
        end
        all_res(frun,:) = options.norm_res;
    end
end
all_res(sum(all_res,2)==0,:) = [];
[~,p]                        = min(sum(all_res,2));
norm_res                     = all_res(p,:);

if strcmp(options.overwrite_data,'on') || ...
        (strcmp(options.overwrite_data,'off') && ~isfield(meta,'normalise'))
    
    fprintf(' running normalize \n')
    if exist('matlabbatch','var')
        clear matlabbatch
    end
    
    % if field maps, normalize EPI from T1
    % -------------------------------------
    if strcmpi(options.norm,'T1norm')
        
        % segment the coregistered T1 (not resliced) and possibly other
        % modality (T2 or FLAIR most likely bit not checking)
        % --------------------------------------------------------
        for channel = 1:size(subject.anat,1)
            if channel == 1
                matlabbatch{1}.spm.spatial.preproc.channel(channel).vols   = {subject.anat{channel}};
            else
                matlabbatch{1}.spm.spatial.preproc.channel(channel).vols   = {others{1}.rfiles{channel-1}};
            end
            matlabbatch{end}.spm.spatial.preproc.channel(channel).biasreg  = 0.001;
            matlabbatch{end}.spm.spatial.preproc.channel(channel).biasfwhm = 60;
            matlabbatch{end}.spm.spatial.preproc.channel(channel).write    = [0 1];
        end
        matlabbatch{end}.spm.spatial.preproc.tissue(1).tpm    = {[spm_root filesep 'tpm' filesep 'TPM.nii,1']};
        matlabbatch{end}.spm.spatial.preproc.tissue(1).ngaus  = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(2).tpm    = {[spm_root filesep 'tpm' filesep 'TPM.nii,2']};
        matlabbatch{end}.spm.spatial.preproc.tissue(2).ngaus  = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(3).tpm    = {[spm_root filesep 'tpm' filesep 'TPM.nii,3']};
        matlabbatch{end}.spm.spatial.preproc.tissue(3).ngaus  = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(4).tpm    = {[spm_root filesep 'tpm' filesep 'TPM.nii,4']};
        matlabbatch{end}.spm.spatial.preproc.tissue(4).ngaus  = 3;
        matlabbatch{end}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(5).tpm    = {[spm_root filesep 'tpm' filesep 'TPM.nii,5']};
        matlabbatch{end}.spm.spatial.preproc.tissue(5).ngaus  = 4;
        matlabbatch{end}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(6).tpm    = {[spm_root filesep 'tpm' filesep 'TPM.nii,6']};
        matlabbatch{end}.spm.spatial.preproc.tissue(6).ngaus  = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.warp.mrf         = 1;
        matlabbatch{end}.spm.spatial.preproc.warp.cleanup     = 1;
        matlabbatch{end}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
        matlabbatch{end}.spm.spatial.preproc.warp.affreg      = 'mni';
        matlabbatch{end}.spm.spatial.preproc.warp.fwhm        = 0;
        matlabbatch{end}.spm.spatial.preproc.warp.samp        = 3;
        if strcmpi(options.norm_inv_save,'on')
            matlabbatch{end}.spm.spatial.preproc.warp.write   = [1 1];
        else
            matlabbatch{end}.spm.spatial.preproc.warp.write   = [0 1];
        end
        
        % normalize EPI using T1 info
        % ----------------------------
        % normalize space-EPI_T1
        matlabbatch{2}.spm.spatial.normalise.write.subj(1).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(1) = {spaceEPI_T1w};
        
        % normalize bias corrected T1 and tissue classes
        matlabbatch{end}.spm.spatial.normalise.write.subj(2).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        
        for channel = 1:size(subject.anat,1)
            matlabbatch{end}.spm.spatial.normalise.write.subj(2).resample(channel) = cfg_dep(['Segment: Bias Corrected (' num2str(channel) ')'], substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
        end
        matlabbatch{end}.spm.spatial.normalise.write.subj(2).resample(channel+1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(2).resample(channel+2) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(2).resample(channel+3) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
        
        % normalize EPI
        matlabbatch{end}.spm.spatial.normalise.write.subj(3).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        for frun = 1:size(realigned_files,1)
            matlabbatch{end}.spm.spatial.normalise.write.subj(3).resample(frun,:) = {realigned_files{frun}}; %#ok<CCAT1>
        end
        matlabbatch{end}.spm.spatial.normalise.write.subj(3).resample(end+1,:)    = {mean_realigned_file}; % adding mean image
        matlabbatch{end}.spm.spatial.normalise.write.woptions.bb                  = [-78 -112 -70 ; 78 76 85];
        matlabbatch{end}.spm.spatial.normalise.write.woptions.vox                 = norm_res;
        matlabbatch{end}.spm.spatial.normalise.write.woptions.interp              = 4;
        matlabbatch{end}.spm.spatial.normalise.write.woptions.prefix              = 'w';
        
    else  % normalize EPI on EPI template (old routine)
        % -------------------------------------------
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source(1)                  = {mean_realigned_file};
        matlabbatch{end}.spm.tools.oldnorm.estwrite.subj.wtsrc                    = '';
        for frun = 1:size(realigned_files,1)
            matlabbatch{end}.spm.tools.oldnorm.estwrite.subj.resample(frun,:)     = {realigned_files{frun}}; %#ok<CCAT1>
        end
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.template = ...
            {[spm_root filesep 'toolbox' filesep 'OldNorm' filesep 'EPI.nii,1']};
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.weight               = '';
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.smosrc               = 8;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.smoref               = 0;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.regtype              = 'mni';
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.cutoff               = 25;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.nits                 = 16;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.reg                  = 1;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.preserve             = 0;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.bb                   = [-78 -112 -70 ; 78 76 85];
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.vox                  = norm_res;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.interp               = 4;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.wrap                 = [0 0 0];
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.prefix               = 'w';
    end
    spm_jobman('run',matlabbatch);
    
    % update json files
    for frun = 1:size(subject.func,1)
        [filepath,filename] = fileparts(subject.func{frun});
        if exist(fullfile(filepath,[filename '.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename '.json']));
        elseif exist(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']));
        elseif exist(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']));
        end
        
        if strcmpi(options.norm,'T1norm')
            meta.segment        = matlabbatch{1}.spm.spatial.preproc.channel.vols;
            meta.normalise.from = 'segment';
        else
            meta.normalise.from = 'EPI template';
        end
        meta.normalise.resolution = norm_res;
        spm_jsonwrite(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']),meta,opts)
    end
    clear matlabbatch;
end

% QC original data
if strcmp(options.QC,'on') && ~isfield(subject,'anat_qa') && ~isfield(anatinfo,'QA')
    [anat_dir,filename] = fileparts(subject.anat{1});
    if contains(filename,'space-IXI549')
        subject.anat_qa = spmup_anatQA(subject.anat{1},...
            fullfile(anat_dir,[filename(1:end-3) 'label-GM_probseg.nii']),...
            fullfile(anat_dir,[filename(1:end-3) 'label-WM_probseg.nii']));
    else
        subject.anat_qa = spmup_anatQA(subject.anat{1},class{1},class{2});
    end
    anatinfo.QA  = subject.anat_qa;
    spm_jsonwrite(fullfile(anat_dir,[filename(1:end-4) '_desc-preprocessed_anat.json']),anatinfo,opts)
end

% ----------
%% Smoothing
% ----------

for frun = 1:size(Normalized_files,1)
    [filepath,filename,ext] = fileparts(Normalized_files{frun});
    stats_ready{frun,1}     = [filepath filesep 's' filename ext];
    
    if strcmp(options.overwrite_data,'on') || ...
            (strcmp(options.overwrite_data,'off') && ~isfield(meta,'smoothingkernel'))
        
        fprintf(' smoothing sesssion/run %g \n',frun);
        V                       = spm_vol(Normalized_files{frun});
        skernel                 = abs(diag(V(1).mat)); clear V
        [filepath,filename,ext] = fileparts(subject.func{frun});
        if exist(fullfile(filepath,[filename '.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename '.json']));
        elseif exist(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']));
        elseif exist(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']),'file')
            meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']));
        end
        
        if strcmp(options.derivatives,'off') || ...
                contains(stats_ready{frun,1},'task-rest','IgnoreCase',true)
            % user left derivatives default but it's resting state
            % (makes no sense to do the other case)
            if isempty(options.skernel)
                options.skernel      = 3;
                meta.smoothingkernel = 3;
            end
            spm_smooth(Normalized_files{frun},stats_ready{frun},skernel(1:3).*options.skernel);
            meta.smoothingkernel = skernel(1:3).*options.skernel;
        else % tiny kernel just to take closest neighbourghs
            spm_smooth(Normalized_files{frun},stats_ready{frun},skernel(1:3));
            meta.smoothingkernel = skernel(1:3);
            if isempty(options.skernel)
                options.skernel = meta.smoothingkernel;
            end
        end
        spm_jsonwrite(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']),meta,opts)
    end
end

% ------------------------------------------------------------------------
%% update the subject structure, add fMRI QC from ready files and clean-up
% ------------------------------------------------------------------------

% sanity check that all images are in the same space.
V_to_check        = stats_ready;
V_to_check{end+1} = Normalized_class{1};
V_to_check{end+1} = Normalized_class{2};
V_to_check{end+1} = Normalized_class{3};
V_to_check{end+1} = NormalizedAnat_file;
if all(cellfun(@(x) exist(x,'file'), V_to_check))
    try
        spm_check_orientations(spm_vol(char(V_to_check)));
    catch
        error('while the spmup_BIDS_preprocess ran and created files, there are oriented differently!! no idea how this happened')
    end
end
clear V_to_check

% since it's all good update file names and the subject structure
% include clean-up as we iterate through files and get QC metrics

% deal with anat
if contains(subject.anat{1},'_T1w.nii')
    [filepath,filename,ext] = fileparts(subject.anat{1});
    partA   = get_fileparts(filename);
    partB   = filename(max(strfind(filename,'_'))+1:end);
    newname = fullfile(filepath,[partA 'space-IXI549_' partB ext]);
    if exist(NormalizedAnat_file,'file')
        movefile(NormalizedAnat_file,newname);
        subject.anat{1} = newname;
    elseif exist(newname,'file')
        subject.anat{1} = newname;
    end
    
    % also rename the json file
    json = dir(fullfile(filepath,[partA '*desc-preprocessed_anat.json']));
    if size(json,1)>1 % repeated stuff, take last one
        [~,index]=max(arrayfun(@(x) x.datenum, json));
        % delete file
        allfiles = [1:size(json,1)];
        allfiles(allfiles == index) = [];
        for d=1:length(allfiles)
           delete(fullfile(json(allfiles(d)).folder,json(allfiles(d)).name));
        end
        %keep right file info
        json = json(index);
    end
    if ~isempty(json)
        newname = fullfile(filepath,[partA 'space-IXI549_desc-preprocessed_anat.json']);
        movefile(fullfile(json.folder,json.name),newname);
    end
else
    error('1st image in subject.anat should be T1w, not sure what to do')
end

labels = {'GM','WM','CSF'};
subject.tissues = Normalized_class';
% indicate in name the prob. segmentation and tissue type
for tissue = 1:3
    if exist(class{tissue},'file')
        [filepath,filename,ext] = fileparts(class{tissue});
        partA = get_fileparts(filename);
        partB = ['label-' labels{tissue} '_probseg' ext];
        movefile(class{tissue},fullfile(filepath,[partA partB]));
    end
    
    [filepath,filename,ext] = fileparts(Normalized_class{tissue});
    partA = get_fileparts(filename);
    partB = ['space-IXI549_label-' labels{tissue} '_probseg' ext];
    if exist(Normalized_class{tissue},'file')
        movefile(Normalized_class{tissue},fullfile(filepath,[partA partB]));
        subject.tissues{tissue} = fullfile(filepath,[partA partB]);
    elseif exist(fullfile(filepath,[partA partB]),'file')
        subject.tissues{tissue} = fullfile(filepath,[partA partB]);
    end
end

% remove bias corrected in T1w space
if strcmpi(options.keep_data,'off')
    msub = dir(fullfile(fileparts(subject.anat{1}),'msub*.nii'));
    if ~isempty(msub)
        delete(fullfile(msub.folder,msub.name))
    end
end

% delete the anat registered to EPI
spaceEPI = dir(fullfile(fileparts(subject.anat{1}),'space-EPI*.nii'));
if ~isempty(spaceEPI)
    delete(fullfile(spaceEPI.folder,spaceEPI.name))
end

rsub = dir(fullfile(fileparts(subject.anat{1}),'rsub*.nii'));
if ~isempty(rsub)
    delete(fullfile(rsub.folder,rsub.name))
end

% rename the normalized anat registered to EPI (the 'true' anat under EPI)
wspaceEPI = dir(fullfile(fileparts(subject.anat{1}),'wspace-EPI*.nii'));
if ~isempty(wspaceEPI)
    partA = get_fileparts(wspaceEPI.name);
    partB = ['space-IXI549_desc-epiresliced' ext];
    movefile(fullfile(wspaceEPI.folder,wspaceEPI.name),fullfile(fileparts(subject.anat{1}),[partA partB]))
end

% delete other coreg images
spaceT1 = dir(fullfile(fileparts(subject.anat{1}),'*space-T1*.nii'));
if ~isempty(spaceT1)
    for f=1:size(spaceT1,1)
        delete(fullfile(spaceT1(f).folder,spaceT1(f).name))
    end
end

% clean-up intermediate fmri files
if strcmpi(options.keep_data,'off')
    for frun = 1:size(subject.func, 1)
        if bold_include(frun)
            [filepath,filename,ext] = fileparts(subject.func{frun});
            if ~isempty(options.removeNvol)
                delete(fullfile(filepath,[filename(1:end-5) '_skip-' num2str(options.removeNvol) '_bold*'])); % volume edited
            end
            if ~strcmp(options.despike,'off')
                delete(fullfile(filepath,[filename(1:end-5) '*_rec-despiked_bold' ext])); % despiked
            end
            st = dir(fullfile(filepath,['st_*' ext]));
            if ~isempty(st)
                for f=1:size(st,1)
                    delete(fullfile(st(f).folder,st(f).name))
                end
            end
            st = dir(fullfile(filepath,'st_*.mat'));
            if ~isempty(st)
                for f=1:size(st,1)
                    delete(fullfile(st(f).folder,st(f).name))
                end
            end
            ur = dir(fullfile(filepath,['ur*' ext]));
            if ~isempty(ur)
                delete(fullfile(ur.folder,ur.name))
            end
            ust = dir(fullfile(filepath,['ust*' ext]));
            if ~isempty(ust)
                delete(fullfile(ust.folder,ust.name))
            end
            wst = dir(fullfile(filepath,['wst*' ext]));
            if ~isempty(wst)
                delete(fullfile(wst.folder,wst.name))
            end
            wur = dir(fullfile(filepath,['wur*' ext]));
            if ~isempty(wur)
                delete(fullfile(wur.folder,wur.name))
            end
            wfmag = dir(fullfile(filepath,['wfmag*' ext]));
            if ~isempty(wfmag)
                delete(fullfile(wfmag.folder,wfmag.name))
            end
            avg = dir(fullfile(filepath,['mean*' ext]));
            if ~isempty(avg)
                delete(fullfile(avg.folder,avg.name))
            end
            wavg = dir(fullfile(filepath,['wmean*' ext]));
            if ~isempty(wavg)
                delete(fullfile(wavg.folder,wavg.name))
            end
        end
    end
end

% rename files and last QC
for frun = 1:size(subject.func,1)
    if exist(stats_ready{frun},'file')
        [filepath,filename,ext] = fileparts(subject.func{frun});
        newname = fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold' ext]);
        movefile(stats_ready{frun},newname);
        subject.func{frun} = newname;
        
        if strcmp(options.QC,'on')
            if exist(fullfile(filepath,[filename '.json']),'file')
                meta = spm_jsonread(fullfile(filepath,[filename '.json']));
            elseif exist(fullfile(filepath,[subject.func{frun}(1:end-9) '_space-IXI549_desc-preprocessed_bold.json']),'file')
                meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']));
            elseif exist(fullfile(filepath,[subject.func{frun}(1:end-9) '_desc-preprocessed_bold.json']),'file')
                meta = spm_jsonread(fullfile(filepath,[subject.func{frun}(1:end-9) '_desc-preprocessed_bold.json']));
            end
            
            if ~isfield(subject.func_qa{frun},'preproc_outliers')
                [~,subject.func_qa{frun}.preproc_outliers] = spmup_volumecorr(newname);
                meta.QA.correlation_outliers = subject.func_qa{frun}.preproc_outliers;
            end
            if ~isfield(subject.func_qa{frun},'preproc_tSNR')
                subject.func_qa{frun}.preproc_tSNR = spmup_temporalSNR(newname,subject.tissues);
                meta.QA.tSNR = subject.func_qa{frun}.preproc_tSNR;
            end
            
            spm_jsonwrite([subject.func{frun}(1:end-4) '.json'],meta,opts)
            oldjson = dir(fullfile(filepath,[filename(1:end-5) '_desc-preprocessed_bold.json']));
            if ~isempty(oldjson)
                delete(fullfile(oldjson.folder,oldjson.name))
            end
        else
            json = dir(fullfile(filepath,[filename(1:end-5) '*desc-preprocessed_bold.json']));
            if ~isempty(json)
                newname = fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']);
                movefile(fullfile(json.folder,json.name),newname);
            end
        end
        
    else % it could be it was done and just need updated
        
        [filepath,filename,ext] = fileparts(subject.func{frun});
        newname = fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold' ext]);
        if exist(newname,'file')
            subject.func{frun} = newname;
            
            if ~isfield(subject,'func_qa')
                subject.func_qa = cell(size(subject.func,1),1);
            end
            
            if strcmp(options.QC,'on') && isempty(subject.func_qa{frun})
                
                if exist(fullfile(filepath,[filename '.json']),'file')
                    meta = spm_jsonread(fullfile(filepath,[filename '.json']));
                elseif exist(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']),'file')
                    meta = spm_jsonread(fullfile(filepath,[filename(1:end-5) '_space-IXI549_desc-preprocessed_bold.json']));
                else
                    meta.QA = [];
                end
                
                if ~isfield(subject.func_qa{frun},'preproc_outliers')
                    if isfield(meta.QA,'correlation_outliers')
                        subject.func_qa{frun}.preproc_outliers  = meta.QA.correlation_outliers;
                    else
                        [~,subject.func_qa{frun}.preproc_outliers] = spmup_volumecorr(newname);
                    end
                end
                
                if ~isfield(subject.func_qa{frun},'preproc_tSNR')
                    if isfield(meta.QA,'tSNR')
                        subject.func_qa{frun}.preproc_tSNR = meta.QA.tSNR;
                    else
                        subject.func_qa{frun}.preproc_tSNR = spmup_temporalSNR(newname,subject.tissues);
                    end
                end
            end
        end
    end
end
end

%% sub-routine

function partA = get_fileparts(filename)
% uses BIDS key-values pairs to reconstruct file name

values  = strfind(filename,'_')-1;
sub     = strfind(filename,'sub-');
sub     = filename(sub:min(values(values>sub))+1);
ses     = strfind(filename,'ses-');
if ~isempty(ses)
    ses  = filename(ses:min(values(values>ses))+1);
end
acq      = strfind(filename,'acq-');
if ~isempty(acq)
    acq  = filename(acq:min(values(values>acq))+1);
end
rec      = strfind(filename,'rec-');
if ~isempty(rec)
    rec  = filename(rec:min(values(values>rec))+1);
end
task     = strfind(filename,'task-');
if ~isempty(task)
    task = filename(task:min(values(values>task))+1);
end
partA = [sub ses acq rec task];
end



