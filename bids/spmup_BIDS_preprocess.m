function [anatQA, fMRIQA, subjects] = spmup_BIDS_preprocess(BIDS, subjects, s, options)

% routine to preprocess BIDS fMRI data following 'options'
%
% FORMAT spmup_BIDS_preprocess(BIDS_dir, BIDS, subjects, s)
%        spmup_BIDS_preprocess(BIDS_dir, BIDS, subjects, s, options)
%
% INPUT - BIDS: the structure returned by spm_BIDS
%       - subjects: a structure containing the fullpath of the unpacked anat,
%           fmap and func files for each subject (see spmup_BIDS_unpack)
%       - s is the subject index to preprocess
%       - options the structure of SPMUP options (see spmup_getoptions)
%
% Example usage: BIDS_dir        = 'F:\WakemanHenson_Faces\fmri';
%                options         = spmup_getoptions(BIDS_dir);
%                options.overwrite_data = 'off';
%                [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options);
%                [anatQA, fMRIQA, subjects] = spmup_BIDS_preprocess(BIDS, subjects, 1, options)
%
% Cyril Pernet & Remi Gau
% --------------------------
%  Copyright (C) SPMUP Team 

%% check a few things
if strcmpi(options.norm,'T1norm') && isempty(options.VDM)
    warning(' Calhoun et al. (2017) showed T1 normalization is not good without FieldMaps \n (%s) add VDM map to options \n switching to ''EPInorm''', ...
        'https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.23737')
    options.norm = 'EPI';
end

% each subject build a job structure around matlabbatch
% this is where each subject is analysed using the various options
spm_root = fileparts(which('spm'));
spm_jobman('initcfg');
opts.indent = ' '; % for spm_jsonwrite
fprintf('\n\nrunning subject %g \n',s)

% ---------------------------------------
% reorient anat file to template
% ---------------------------------------

RMexist = 0;
[target_dir,filename] = fileparts(subjects{s}.anat);
if exist(fullfile(target_dir,[filename '_preprocessed.json']),'file')
    anatinfo = spm_jsonread(fullfile(target_dir,[filename '_preprocessed.json']));
    if isfield(anatinfo,'reorientation_matrix')
        RMexist = 1;
    end
end
    
if strcmp(options.overwrite_data,'on') || ...
        (strcmp(options.overwrite_data,'off') && ~RMexist)
    
    Data = [subjects{s}.anat; subjects{s}.func];
    if isfield(subjects{s},'fieldmap')
        Data = [Data ; subjects{s}.fieldmap];
    end
    if ~isempty(options.VDM)
        Data = [Data ; options.VDM{s}]; % this assumes only 1 field map though, will need to be changed
    end
    RM = spmup_auto_reorient(Data,1);
    disp(' Data reoriented');
   
    % update info about all those files
    for f=1:size(Data,1)
       [filepath,filename] = fileparts(Data{f});
       if f==1
           meta.reoriented = '(0,0,0) set with spmup_auto_reorient';
           meta.reorientation_matrix = RM;
       else
           meta.reoriented = '(0,0,0) set from anat with spmup_auto_reorient';
       end
       spm_jsonwrite(fullfile(filepath,[filename '_preprocessed.json']),meta,opts)
       clear meta
    end
end

% ---------------------------------------
% Despiking and slice timing for each run
% ----------------------------------------
bold_include = [];
included_idx = 0;
for frun = 1:size(subjects{s}.func, 1) % each run
    
    filesin = subjects{s}.func{frun};
    [filepath,filename,ext] = fileparts(filesin);
    if exist(fullfile(filepath,[filename '_preprocessed.json']),'file')
        meta = spm_jsonread(fullfile(filepath,[filename '_preprocessed.json'])); %#ok<NASGU>
    end
    metadata.run = frun;
    metadata.TaskName = subjects{s}.func_metadata{frun}.TaskName;
    
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
        hdr          = spm_vol(subjects{s}.func{frun});
        epi_res      = diag(hdr(1).mat);
        epi_res(end) = [];
        
        if ~isfield(subjects{s}.func_metadata{frun}, 'SliceTiming')
            error('Slice Timing Information missing')
        else
            SliceTiming = subjects{s}.func_metadata{frun}.SliceTiming;
        end
        
        if ~isfield(subjects{s}.func_metadata{frun}, 'RepetitionTime')
            error('RepetitionTime Information missing')
        else
            RepetitionTime = subjects{s}.func_metadata{frun}.RepetitionTime;
        end
        
        % -----------------------------------------
        % remove first volumes
        % -----------------------------------------
        
        if options.removeNvol ~=0
            if strcmp(options.overwrite_data,'on')...
                    || ((strcmp(options.overwrite_data,'off') && ~isfield(meta,'NumberOfVolumesDiscarded')))
                
                % remove from the 4D the files we don't want and proceed
                fprintf('\nadjusting 4D file sizes run %g \n',frun)
                
                three_dim_files         = spm_file_split(filesin);
                V                       = three_dim_files; 
                V(1:options.removeNvol) = [];
                spm_file_merge(V,filesin);
                spm_unlink(three_dim_files.fname)
                
                % write a json file containing the details of what volumes were
                % removed (see BIDS derivatives specs)
                meta.NumberOfVolumesDiscarded = options.removeNvol;
                spm_jsonwrite(fullfile(filepath,[filename '_preprocessed.json']),meta,opts)
            end
        end
        
        % -----------------------------------------
        % despiking using adaptive median filter
        % -----------------------------------------
        
        if strcmp(options.despike,'before')
            if strcmp(options.overwrite_data,'on') || ...
                    (strcmp(options.overwrite_data,'off') && ~isfield(meta,'Despiked'))
                flags = struct('auto_mask','on', 'method','median', 'window', [],'skip',0, 'savelog', 'off');
                [Vin,~,loginfo] = spmup_despike(fullfile(filepath,[filename,ext]),[],flags);
                filesin = Vin.fname; clear Vin;

                meta.Despiked = 'before slice timing using spmup_despike';
                meta.despike.param = 'median filter based on the autocorrelation function';
                meta.despike.RMS   = loginfo.RMS;
                spm_jsonwrite(fullfile(filepath,[filename '_preprocessed.json']),meta,opts)
            else
                if exist(fullfile(filepath, ['despiked_' filename ext]),'file')
                    filesin = fullfile(filepath, ['despiked_' filename ext]);
                else
                    filesin = fullfile(filepath, [filename ext]);
                end
            end
            [filepath,filename,ext] = fileparts(filesin);
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
            
            st_files{included_idx,1} = [filepath filesep 'st_' filename ext]; %#ok<*AGROW>
            
            if strcmp(options.overwrite_data, 'on') || (strcmp(options.overwrite_data, 'off') ...
                    && ~exist(st_files{included_idx,1}, 'file'))
                fprintf('\n\nstarting slice timing correction run %g subject %g \n',frun,s)
                spm_slice_timing(filesin, sliceorder, refslice, timing, 'st_');
            end
        else
            st_files{included_idx,1} = fullfile(filepath, [filename ext]);
        end
    end
end % end processing per run


% ----------------------
% Field Map - compute VDM
% ----------------------

if strcmp(options.ignore_fieldmaps, 'off') && isfield(subjects{s}, 'fieldmap')
    
    for ifmap = 1:numel(subjects{s}.fieldmap)
        
        % list which fieldmap goes to which bold run and vice versa
        [~, filename] = spm_fileparts(subjects{s}.fieldmap(ifmap).metadata.IntendedFor);
        which_func_file = contains(subjects{s}.func, filename);
        subjects{s}.fieldmap(ifmap).metadata.which_func = find(which_func_file); % to know which func files needs this fmap
        subjects{s}.which_fmap(find(which_func_file)) = ifmap; % to know which fmap is needed by which func files
        
        % only continue if this fmap is needed for any of the bold run for
        % this task / acq / rec
        if any(bold_include(which_func_file))
            
            % find bold reference run to which the fielmap will be coregistered
            % the first of that task / acq / rec
            bold_ref = subjects{s}.func{find(bold_include, 1)};
            [filepath,filename,ext]=fileparts(bold_ref);
            
            % computes average of that run
            if strcmp(options.despike,'on')
                avg = [filepath filesep 'spmup_mean_st_despiked_' filename ext];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                        && ~exist(avg,'file'))
                    spmup_bascis([filepath filesep 'st_despiked_' filename ext], 'mean'); % use the mean despiked slice timed EPI for QC
                end
            else
                avg = [filepath filesep 'spmup_mean_st_' filename ext];
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                        && ~exist(avg,'file'))
                    spmup_bascis([filepath filesep 'st_' filename ext],'mean'); % use the mean slice timed EPI for QC
                end
            end
            
            % VDM creation: output images
            % check which types of fieldmaps we are dealing with
            switch subjects{s}.fieldmap(ifmap).type
                case 'phasediff'
                    % two magnitudes (use only 1) and 1 phase difference image
                    [filepath,name,ext]=fileparts(subjects{s}.fieldmap(ifmap).phasediff);
                    vdm{ifmap,1} = [filepath filesep 'vdm5_sc' name ext];
                case 'phase12'
                    % two magnitude images and 2 phase images
                    [filepath,name,ext]=fileparts(subjects{s}.fieldmap(ifmap).phase1);
                    vdm{ifmap,1} = [filepath filesep 'vdm5_sc' name ext];
                case 'fieldmap'
                    warning('Fieldmap type of fielmaps not implemented')
                case 'epi'
                    warning('EPI type of fielmaps not implemented')
                otherwise
                    warning('%s is an unsupported type of fieldmap', subjects{s}.fieldmap(ifmap).type)
            end
            
            
            if strcmp(subjects{s}.fieldmap(ifmap).type, 'phasediff') || ...
                    strcmp(subjects{s}.fieldmap(ifmap).type, 'phase12')
                
                if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                        && ~exist(vdm{ifmap},'file'))
                    
                    % coregister fmap to bold reference
                    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {avg};
                    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep = [16 8 4 2 1];
                    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol = ...
                        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
                    
                    if strcmp(subjects{s}.fieldmap(ifmap).type, 'phasediff')
                        % two magnitudes (use only 1) and 1 phase difference image
                        matlabbatch{end}.spm.spatial.coreg.estimate.source = ...
                            {subjects{s}.fieldmap(ifmap).mag1};
                        matlabbatch{end}.spm.spatial.coreg.estimate.other = ...
                            {subjects{s}.fieldmap(ifmap).phasediff};
                        
                    elseif strcmp(subjects{s}.fieldmap(ifmap).type, 'phase12')
                        % two magnitide images and 2 phase images
                        matlabbatch{end}.spm.spatial.coreg.estimate.source = ...
                            {subjects{s}.fieldmap(ifmap).mag1};
                        matlabbatch{end}.spm.spatial.coreg.estimate.other = {...
                            subjects{s}.fieldmap(ifmap).mag2;...
                            subjects{s}.fieldmap(ifmap).phase1;...
                            subjects{s}.fieldmap(ifmap).phase2};
                    end
                    
                    fprintf('\ncoregistering fmap %g subject %g to its reference run ',ifmap,s)
                    spm_jobman('run',matlabbatch); clear matlabbatch;
                    
                    
                    % VDM creation: input images
                    if strcmp(subjects{s}.fieldmap(ifmap).type, 'phasediff')
                        % two magnitudes (use only 1) and 1 phase difference image
                        matlabbatch = get_FM_workflow('phasediff');
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = ...
                            {subjects{s}.fieldmap(ifmap).phasediff};
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = ...
                            {subjects{s}.fieldmap(ifmap).mag1};
                    elseif strcmp(subjects{s}.fieldmap(ifmap).type, 'phase12')
                        % two magnitide images and 2 phase images
                        matlabbatch = get_FM_workflow('phase&mag');
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase = ...
                            {subjects{s}.fieldmap(ifmap).phase1};
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag   = ...
                            {subjects{s}.fieldmap(ifmap).mag1};
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase  = ...
                            {subjects{s}.fieldmap(ifmap).phase2};
                        matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag    = ...
                            {subjects{s}.fieldmap(ifmap).mag2};
                    end
                    
                    % update parameters
                    
                    % fieldmap metadata
                    fmap_metadata = subjects{s}.fieldmap(ifmap).metadata;
                    echotimes =  1000.*[fmap_metadata.EchoTime1 ...
                        fmap_metadata.EchoTime2]; % change to ms
                    matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = echotimes;
                    
                    % func run metadata: if the fmap is applied to several
                    % runs we take the metadata of the first run it must
                    % applied to
                    func_metadata = subjects{s}.func_metadata{which_func_file};
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
                        if sum(findstr(fmap_metadata.PulseSequenceType,'EPI')) ~= 0
                            matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
                        else
                            matlabbatch{end}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
                        end
                    else
                        disp('using default sequence! assuming non-EPI acquisition')
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
                    
                    fprintf('\ncomputing voxel displacement map %g subject %g',ifmap,s)
                    spm_jobman('run',matlabbatch); clear matlabbatch;
                    
                end
                
            end
        end
    end
end


% -----------------------
% Realignment across runs
% ------------------------

if strcmp(options.realign_unwarp, 'off')
    
    % file out of realign will be
    [filepath,filename,ext] = fileparts(st_files{1});
    mean_realigned_file = [filepath filesep 'mean' filename ext];
    for frun = 1:size(st_files,1)
        realigned_files{frun,1} = st_files{frun}; % because we don't reslice, simple encode the linear transform in the header
        [filepath,filename] = fileparts(st_files{frun});
        multi_reg{frun,1} = [filepath filesep 'rp_' filename '.txt'];
    end
    
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~exist(mean_realigned_file,'file'))
        
        fprintf('subject %g: starting realignment \n',s)
        
        for frun = 1:size(st_files,1)
            matlabbatch{1}.spm.spatial.realign.estwrite.data{frun} = {st_files{frun}};
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
        
        spm_jobman('run',matlabbatch); clear matlabbatch;
        
    end
else
    
    %which fieldmap for each run of this task / acq /rec
    which_fmap = subjects{s}.which_fmap(find(bold_include));
    
    % file out of realign will be
    [filepath,filename,ext] = fileparts(st_files{1});
    mean_realigned_file = [filepath filesep 'meanur' filename ext];
    for frun = 1:size(st_files,1)
        [filepath,filename,ext] = fileparts(st_files{frun});
        realigned_files{frun,1} = [filepath filesep 'ur' filename ext]; % because we have the reslice here (not linear)
        multi_reg{frun,1} = [filepath filesep 'rp_' filename '.txt'];
    end
    
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~exist(mean_realigned_file,'file'))
        
        fprintf('subject %g: starting realignment and unwarping \n',s)
        
        for frun = 1:size(st_files,1)
            matlabbatch{1}.spm.spatial.realignunwarp.data(frun).scans = {st_files{frun}};
            if strcmp(options.ignore_fieldmaps, 'off') && isfield(subjects{s}, 'fieldmap')
                matlabbatch{end}.spm.spatial.realignunwarp.data(frun).pmscan = {vdm{which_fmap(frun)}};
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
        
        spm_jobman('run',matlabbatch); clear matlabbatch;
        
        
        % cleanup
        if strcmpi(options.keep_data,'off')
            delete(st_files);
        end
        
    end
    
end


% ---------------
% Coregistration / Segmentation
% ---------------

[filepath,filename,ext] = fileparts(subjects{s}.anat);
EPI_class{1} = [filepath filesep 'c1r' filename ext];
EPI_class{2} = [filepath filesep 'c2r' filename ext];
EPI_class{3} = [filepath filesep 'c3r' filename ext];
if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
        && ~exist(EPI_class{1},'file'))
    
    fprintf('\n\nsubject %g: coregister, segment \n',s)
    
    if exist('matlabbatch','var')
        clear matlabbatch
    end
    
    % coregister anatomical to mean EPI
    % ---------------------------------
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {mean_realigned_file};
    matlabbatch{end}.spm.spatial.coreg.estimate.source = {subjects{s}.anat};
    matlabbatch{end}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep = [16 8 4 2 1];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol = ...
        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    % reslice anatomical to mean EPI
    matlabbatch{2}.spm.spatial.coreg.write.ref = {mean_realigned_file};
    matlabbatch{end}.spm.spatial.coreg.write.source = {subjects{s}.anat};
    matlabbatch{end}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{end}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{end}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{end}.spm.spatial.coreg.write.roptions.prefix = 'r';
    
    
    % run the segmentation on the resliced T1 to get tissue classes
    % in the same space as the EPI data before mormalization
    % -------------------------------------------------------
    matlabbatch{3}.spm.spatial.preproc.channel.vols(1) = ...
        cfg_dep('Coregister: Reslice: Resliced Images', ...
        substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','rfiles'));
    matlabbatch{end}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{end}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{end}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(1).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,1']};
    matlabbatch{end}.spm.spatial.preproc.tissue(1).ngaus = 2;
    matlabbatch{end}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(2).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,2']};
    matlabbatch{end}.spm.spatial.preproc.tissue(2).ngaus = 2;
    matlabbatch{end}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(3).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,3']};
    matlabbatch{end}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{end}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(4).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,4']};
    matlabbatch{end}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{end}.spm.spatial.preproc.tissue(4).native = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(5).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,5']};
    matlabbatch{end}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{end}.spm.spatial.preproc.tissue(5).native = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(6).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,6']};
    matlabbatch{end}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{end}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{end}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{end}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{end}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{end}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{end}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{end}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{end}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{end}.spm.spatial.preproc.warp.write = [0 0];
    
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
    
end


% ---------------
% Normalization
% ---------------

[filepath,filename,ext] = fileparts(subjects{s}.anat);
% Normalization_file = [filepath filesep 'y_' filename ext];
NormalizedAnat_file  = [filepath filesep 'wm' filename ext];
% NormalizedReslicedAnat_file  = [filepath filesep 'wr' filename ext];
class{1} = [filepath filesep 'c1' filename ext];
class{2} = [filepath filesep 'c2' filename ext];
class{3} = [filepath filesep 'c3' filename ext];
Normalized_class{1} = [filepath filesep 'wc1' filename ext];
Normalized_class{2} = [filepath filesep 'wc2' filename ext];
Normalized_class{3} = [filepath filesep 'wc3' filename ext];
for frun = 1:size(realigned_files,1)
    [filepath,filename,ext] = fileparts(realigned_files{frun});
    Normalized_files{frun,1} = [filepath filesep 'w' filename ext];
end

if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
        && ~exist(Normalized_files{end},'file'))
    
    fprintf('\n\nsubject %g: normalize \n',s)
    
    if exist('matlabbatch','var')
        clear matlabbatch
    end
    
    % if field maps, normalize EPI from T1
    % -------------------------------------
    if strcmpi(options.norm,'T1norm')
        
        % segment the coregistered T1 (not resliced)
        % -------------------------------------------
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {subjects{s}.anat};
        matlabbatch{end}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{end}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{end}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{end}.spm.spatial.preproc.tissue(1).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,1']};
        matlabbatch{end}.spm.spatial.preproc.tissue(1).ngaus = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(2).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,2']};
        matlabbatch{end}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(3).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,3']};
        matlabbatch{end}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(4).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,4']};
        matlabbatch{end}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{end}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(5).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,5']};
        matlabbatch{end}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{end}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(6).tpm = {[spm_root filesep 'tpm' filesep 'TPM.nii,6']};
        matlabbatch{end}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{end}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{end}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{end}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{end}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{end}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{end}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{end}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{end}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{end}.spm.spatial.preproc.warp.write = [1 1];
        
        % normalize EPI using T1 info
        % ----------------------------
        % normalize T1
        matlabbatch{2}.spm.spatial.normalise.write.subj(1).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(1) = ...
            cfg_dep('Segment: Bias Corrected (1)', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(2) = ...
            cfg_dep('Segment: c1 Images', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(3) = ...
            cfg_dep('Segment: c2 Images', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
        matlabbatch{end}.spm.spatial.normalise.write.subj(1).resample(4) = ...
            cfg_dep('Segment: c3 Images', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
        
        % normalize EPI
        matlabbatch{end}.spm.spatial.normalise.write.subj(2).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        for frun = 1:size(realigned_files,1)
            matlabbatch{end}.spm.spatial.normalise.write.subj(2).resample(frun,:) = {realigned_files{frun}};
        end
        matlabbatch{end}.spm.spatial.normalise.write.subj(2).resample(end+1,:) = {mean_realigned_file}; % adding mean image
        
        matlabbatch{end}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70 ; 78 76 85];
        matlabbatch{end}.spm.spatial.normalise.write.woptions.vox = epi_res;
        matlabbatch{end}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{end}.spm.spatial.normalise.write.woptions.prefix = 'w';
        
    else  % normalize EPI on EPI template (old routine)
        % -------------------------------------------
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source(1) = {mean_realigned_file};
        matlabbatch{end}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
        for frun = 1:size(realigned_files,1)
            matlabbatch{end}.spm.tools.oldnorm.estwrite.subj.resample(frun,:) = {realigned_files{frun}};
        end
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.template = ...
            {[spm_root filesep 'toolbox' filesep 'OldNorm' filesep 'EPI.nii,1']};
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -70 ; 78 76 85];
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.vox = epi_res;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.interp = 4;
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{end}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
        
    end
    
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
    
end


% --------------
% Smoothing
% --------------

if strcmp(options.derivatives,'off') % otherwise do it after stats
    
    for frun = 1:size(subjects{s}.func,1)
        [filepath,filename,ext] = fileparts(Normalized_files{frun});
        stats_ready{frun,1} = [filepath filesep 's' filename ext];
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                && ~exist(stats_ready{frun},'file'))
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
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~exist([fileparts(NormalizedAnat_file) filesep 'anatQA.mat'],'file'))
        
        fprintf('subject %g: Anatomical Quality control \n',s)
        % sanity check that all images are in the same space.
        V_to_check = Normalized_class';
        V_to_check{end+1} = NormalizedAnat_file;
        spm_check_orientations(spm_vol(char(V_to_check)));
        
        % Basic QA for anatomical data is to get SNR, CNR, FBER and Entropy
        % This is useful to check coregistration and normalization worked fine
        tmp = spmup_anatQA(NormalizedAnat_file,Normalized_class{1},Normalized_class{2});
        save([fileparts(NormalizedAnat_file) filesep 'anatQA.mat'],'tmp');
        
    else
        load([fileparts(NormalizedAnat_file) filesep 'anatQA.mat'],'tmp');
    end
    
    anatQA = tmp; clear tmp
    
    if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
            && ~exist([fileparts(Normalized_files{end}) filesep 'fMRIQA.mat'],'file') )
        
        fprintf('subject %g: fMRI Quality control \n',s)
        % For functional data, QA is consists in getting temporal SNR and then
        % check for motion - here we also compute additional regressors to
        % account for motion
        
        davg = spmup_comp_dist2surf(subjects{s}.anat);
        
        if strcmpi(options.scrubbing,'on')
            flags = struct('motion_parameters','on','globals','on','volume_distance','off','movie','off', ...
                'AC', [], 'average','on', 'T1', 'on');
        else
            flags = struct('motion_parameters','off','globals','off','volume_distance','on','movie','off', ...
                'AC', [], 'average','on', 'T1', 'on');
        end
        
        for frun = 1:size(stats_ready,1)
            
            % sanity check that all images are in the same space.
            V_to_check = Normalized_class';
            V_to_check{end+1} = stats_ready{frun};
            spm_check_orientations(spm_vol(char(V_to_check)));
            
            %             fMRIQA.tSNR(1, frun) = spmup_temporalSNR(Normalized_files{frun}, Normalized_class, 1);
            fMRIQA.tSNR(1, frun) = spmup_temporalSNR(Normalized_files{frun}, Normalized_class, 0);
            
            tmp = spmup_first_level_qa(NormalizedAnat_file, cell2mat(stats_ready(frun)), flags);
            fMRIQA.meanFD(1,frun) = mean(spmup_FD(cell2mat(tmp), davg));
            clear tmp
            
            QA.tSNR = fMRIQA.tSNR(1,frun);
            QA.meanFD = fMRIQA.meanFD(1,frun);
            
            save([fileparts(Normalized_files{frun}) filesep 'fMRIQA.mat'],'QA');
            clear QA
            
        end
        
    else
        load([fileparts(Normalized_files{frun}) filesep 'fMRIQA.mat'],'QA');
        
        fMRIQA.tSNR(1, frun) = QA.tSNR;  %#ok<NODEF>
        fMRIQA.meanFD(1,frun) = QA.meanFD;
        clear QA
    end
else
    anatQA = [];
    fMRIQA.tSNR = [];
    fMRIQA.meanFD = [];
end

if strcmp(options.carpet_plot,'on')
    % create carpet plots
    for frun = 1:size(subjects{s}.func, 1)
        if bold_include(frun)
            if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                    && ~exist(fullfile(fileparts(realigned_files{frun}), 'voxplot.fig'), 'file'))
                fprintf('subject %g: fMRI Quality control: carpet plot \n',s)
                P = subjects{s}.func{frun};
                c1 = EPI_class{1};
                c2 = EPI_class{2};
                c3 = EPI_class{3};
                spmup_timeseriesplot(P, c1, c2, c3, 'motion','on','nuisances','on','correlation','on');
            end
        end
    end
end



% Clean up
if strcmpi(options.keep_data,'off')
    
    for frun = 1:size(subjects{s}.func, 1)
        if bold_include(frun) 
            if strcmp(options.despike,'on')
                delete(subjects{s}.func{frun}); % original
            end
            if strcmp(options.despike,'on') && exist('SliceTiming', 'var')
                [filepath,filename,ext]=fileparts(subjects{s}.func{frun});
                delete(fullfile(filepath, ['despiked_' filename ext])); % despiked
            end
        end
    end
    
    for frun = size(realigned_files,1)
        delete(realigned_files{frun});
        if strcmp(options.derivatives,'off')
            delete(Normalized_files{frun});
        end
    end
end

end

%% field map sub-routine
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

