function subject = spmup_BIDS_1rstlevel(subject, options)

% 1st level analysis of BIDS fMRI data
%
% FORMAT subject = spmup_BIDS_1rstlevel(subject, options)
%
% INPUT - subject: a structure containing the fullpath of the unpacked anat,
%           fmap and func files for a given subject (see spmup_BIDS_unpack)
%       - options the structure of SPMUP options (see spmup_getoptions)
%
% OUTPUT subject is the same structure updated
%        - with the preprocessed data information
%        - with QC metrics
%
% Example usage from start to finish: BIDS_dir        = 'F:\WakemanHenson_Faces\fmri';
%                                     options         = spmup_getoptions(BIDS_dir);
%                                     [BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options);
%                                     subject         = spmup_BIDS_preprocess(subjects{n}, options)
%                                     subject         = spmup_BIDS_1rstlevel(subjects{n}, options)
%
%         usage after preprocesing is done, i.e. just run stats modelling
% 
%
% Cyril Pernet 
% --------------------------
%  Copyright (C) SPMUP Team 

opts.indent = ' '; % for spm_jsonwrite
spm_jobman('initcfg');

%% get the 1st level motion regressors updated

% if possible get the right distance between (0,0,0) and the surface
[filepath,filename,ext] = fileparts(subject.anat{1});
if exist(fullfile(filepath,[filename(1:end-4) '_desc-preprocessed_anat.json']),'file')
    anatinfo = spm_jsonread(fullfile(filepath,[filename(1:end-4) '_desc-preprocessed_anat.json']));
    if isfield(anatinfo,'distance2surface')
        davg = anatinfo.distance2surface;
    end
end
    
if ~exist('davg','var')
    c1 = dir(fullfile(filepath,[filename(1:end-4) '_label-GM_probseg' ext]));
    if isempty(c1)
        c1 = dir(fullfile(filepath,['c1' filename ext]));
    end
    c2 = dir(fullfile(filepath,[filename(1:end-4) '_label-WM_probseg' ext]));
    if isempty(c1)
        c1 = dir(fullfile(filepath,['c2' filename ext]));
    end
    
    if ~isempty(c1) && ~isempty(c2)
        c1   = fullfile(c1.folder,c1.name);
        c2   = fullfile(c2.folder,c2.name);
        davg = spmup_comp_dist2surf(subject.anat{1},c1,c2);
        anatinfo.distance2surface = davg;
    else
        davg = 50;
    end
end

if davg ~= 50
    spm_jsonwrite(fullfile(filepath,[filename(1:end-4) '_desc-preprocessed_anat.json']),anatinfo,opts)
end

% now compute
for frun = max(size(subject.func)):-1:1  % max of cell array so it doesn't matter transposed or not
    if strcmpi(options.scrubbing,'on')
        QAjobs{frun} = spmup_first_level_qa(subject.func{frun}, 'Radius',davg, ...
            'Movie','off','Voltera',options.motionexp,...
            'Framewise displacement','on','Globals','on');
    else
        QAjobs{frun} = spmup_first_level_qa(subject.func{frun}, 'Radius',davg, ...
            'Movie','off','Voltera',options.motionexp,...
            'Framewise displacement','off','Globals','off');
    end
    
    if contains(subject.func{frun},'task-rest','IgnoreCase',true)
        design = load(QAjobs{frun}.design);
        json   = jsondecode(fileread([QAjobs{frun}.design(1:end-4) '.json']));
        
        if strcmpi(options.nuisance,'WMCSF')
            QAjobs{frun}.nuisance = spmup_nuisance(subject.func{frun},subject.tissues{2},subject.tissues{3});
            design = [design QAjobs{frun}.nuisance.WM QAjobs{frun}.nuisance.CSF]; %#ok<AGROW>
            json.Columns{end+1} = 'WM';
            json.Columns{end+1} = 'CSF';
        elseif strcmpi(options.nuisance,'compcor')
            brainmask             = spmup_auto_mask(spm_vol(subject.func{frun}));
            [anoise,tnoise]       = spmup_rcompcor(subject.func{frun},brainmask,subject.tissues{2},...
                subject.tissues{3},subject.tissues{1});
            QAjobs{frun}.nuisance.anoise = anoise;
            QAjobs{frun}.nuisance.tnoise = tnoise;
            design = [design anoise{3} tnoise]; %#ok<AGROW>
            json.Columns{end+size(QAjobs{frun}.nuisance.anoise,2)} = 'anoise';
            json.Columns{end+size(QAjobs{frun}.nuisance.tnoise,2)} = 'tnoise';
        else
           error('unknown nuisance to clean resting state data') 
        end
        save(QAjobs{frun}.design,'design','-ascii'); clear design
        spm_jsonwrite([QAjobs{frun}.design(1:end-4) '.json'],json,opts); clear json
    end
    
    if ~isfield(subject.func_qa{frun},'FramewiseDisplacement') && isfield(QAjobs{frun},'FD')
        subject.func_qa{frun}.FD = QAjobs{frun}.FD;
    end
    if ~isfield(subject.func_qa{frun},'Globals') && isfield(QAjobs{frun},'glo')
        subject.func_qa{frun}.Globals = QAjobs{frun}.glo;
    end
    
    if isfield(QAjobs{frun},'nuisance')
        if ~isfield(subject.func_qa{frun},'nuisance')
            subject.func_qa{frun}.nuisance = QAjobs{frun}.nuisance;
        else
            if all(~isfield(subject.func_qa{frun}.nuisance,{'WM','CSF'})) && ...
                all(isfield(QAjobs{frun}.nuisance,{'WM','CSF'})) 
                subject.func_qa{frun}.nuisance = QAjobs{frun}.nuisance;
            elseif ~isfield(subject.func_qa{frun}.nuisance,'anoise') && ...
                    isfield(QAjobs{frun}.nuisance,'anoise')
                subject.func_qa{frun}.nuisance.anoise = QAjobs{frun}.nuisance.anoise;
                subject.func_qa{frun}.nuisance.tnoise = QAjobs{frun}.nuisance.tnoise;
            end
        end
    end
end

% --------------------
% 1st level modelling
% --------------------
filepath    = fileparts(fileparts(subject.func{1}));
Statspath   = [filepath filesep 'spmup_stats'];
SPMmat_file = [Statspath filesep 'SPM.mat'];
if strcmp(options.overwrite_data,'on') || ...
        (strcmp(options.overwrite_data,'off') && ~exist(SPMmat_file,'file'))
    
    clear matlabbatch; 
    if ~exist(Statspath,'dir')
        mkdir(Statspath)
        check_mat = dir(fullfile(Statspath,'SPM.mat'));
        if ~isempty(check_mat)
           delete(fullfile(check_mat.folder,check_mat.name)) 
        end
    end
    N                                                 = length(subject.func_metadata{1}.SliceTiming);
    matlabbatch{1}.spm.stats.fmri_spec.dir            = {Statspath};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units   = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = subject.func_metadata{1}.RepetitionTime;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = N;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(N/ 2);
    
    % sessions, onsets and durations
    for frun = max(size(subject.func)):-1:1 
        if ~contains(subject.func{frun},'task-rest','IgnoreCase',true)
            
            events = readtable(subject.event{frun}, ...
                'FileType', 'text', ...
                'Delimiter', '\t', ...
                'TreatAsEmpty', {'N/A','n/a'});
            
            if contains(options.conditions,events.Properties.VariableNames)
                cond = unique(events.(options.conditions));
            else % note the default options.conditions is trial_type
                warning('%s not find in events.tsv turning back to defaults trial_type',options.conditions)
                if contains('trial_type',events.Properties.VariableNames)
                    cond = unique(events.trial_type);
                    options.conditions = 'trial_type';
                    if size(cond,1)==1 && contains('event_type',events.Properties.VariableNames)
                        warning('number of unique condition found = 1, using stim_type instead') 
                        cond = unique(events.event_type);
                        options.conditions = 'event_type';
                    elseif size(cond,1)==1 && contains('stim_type',events.Properties.VariableNames)
                        warning('number of unique condition found = 1, using stim_type instead') 
                        cond = unique(events.stim_type);
                        options.conditions = 'stim_type';
                    end
                elseif contains('event_type',events.Properties.VariableNames)
                    warning('trial_type not find in events.tsv using event_type instead') 
                    cond = unique(events.event_type);
                    options.conditions = 'event_type';
                elseif contains('stim_type',events.Properties.VariableNames)
                    warning('trial_type not find in events.tsv using stim_type instead') 
                    cond = unique(events.stim_type);
                    options.conditions = 'stim_type';
                end
            end
            % make sure no empty
            rm_cond = zeros(1,length(cond));
            for n=length(cond):-1:1
                if strcmpi(cond{n},'n/a') || isempty(cond{n})
                    rm_cond(n) = 1;
                end
            end
            cond(find(rm_cond)) = []; %#ok<FNDSB>
            N_cond = size(cond,1);
            if N_cond == 1
               warning('number of unique condition found = 1') 
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans = {subject.func{frun}};  %#ok<CCAT1>
            for C = 1:N_cond
                trial_index = cellfun(@(x) strcmp(cond{C},x), events.(options.conditions));
                onsets      = events.onset(trial_index);
                durations   = events.duration(trial_index);
                matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).name     = cond{C};
                matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).onset    = onsets(~isnan(onsets));
                matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).duration = durations(~isnan(durations));
                matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).tmod     = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).orth =    1;
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).multi         = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).regress       = struct('name', {}, 'val', {});
            if strcmpi(options.scrubbing,'off')
                matlabbatch{1}.spm.stats.fmri_spec.sess(frun).multi_reg = {subject.motionfile{frun}}; %#ok<CCAT1>
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(frun).multi_reg = {QAjobs{frun}.design};
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).hpf           = 128;
        else % resting state
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans         = {subject.func{frun}}; %#ok<CCAT1>
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond          = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).multi         = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).regress       = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).multi_reg     = {QAjobs{frun}.design};
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).hpf           = 1000;
        end
    end
    
    % convolution
    matlabbatch{1}.spm.stats.fmri_spec.fact                = struct('name', {}, 'levels', {});
    if options.derivatives == 1 || strcmp(options.derivatives, '1')
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
    elseif options.derivatives == 2 || strcmp(options.derivatives, '2')
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
    else
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    end
    
    % filters, etc .. 
    matlabbatch{1}.spm.stats.fmri_spec.volt     = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global   = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh  = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask     = {''};
    if all(contains(matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans,'rest'))
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'None';
    else
        if subject.func_metadata{frun}.RepetitionTime >=1
            matlabbatch{1}.spm.stats.fmri_spec.cvi  = 'AR(1)';
        else
            matlabbatch{1}.spm.stats.fmri_spec.cvi  = 'FAST';
        end
    end
    
    % estimate
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', ...
        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','spmmat'));
    if strcmp(options.carpet_plot,'on') || ...
            all(contains(matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans,'rest'))
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
    else
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    end
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    save([Statspath filesep 'stats_batch.mat'],'matlabbatch');
    if exist(fullfile(Statspath,'SPM.mat'),'file')
        delete(fullfile(Statspath,'SPM.mat'))
    end
    
    if all(contains(matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans,'rest'))
        residuals = spm_jobman('run',matlabbatch);
        if size(matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans,1) && ...
                ~iscell(matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans)
            [filepath,filename,ext]=fileparts(matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans);
            VfMRI = spm_file_merge(residuals{2}.res,fullfile(filepath,[filename(1:end-5) '-GLMdenoised_bold' ext]));
        else
            [~,filename,ext]=fileparts(matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans{1});
            VfMRI = spm_file_merge(residuals{2}.res,fullfile(fileparts(Statspath),[filename(1:end-5) '-GLMdenoised_bold' ext]));
        end
        subject.func{1} = VfMRI(1).fname;
        spmup_timeseriesplot(VfMRI(1).fname, ...
            subject.tissues{1}, subject.tissues{2}, subject.tissues{3}, ...
            'motion','on','nuisances','on','correlation','off',...
            'clean','off','figure','save');
        
        if strcmpi(options.keep_data,'off')
            rmdir(Statspath,'s');
        else
            if size(matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans,1)
                movefile(Statspath,fullfile(filepath,'GLMdenoise'));
            else
                movefile(Statspath,fullfile(fileparts(Statspath),'GLMdenoise'));
            end
        end
        clear matlabbatch; resting_state = 'done'; %#ok<NASGU>
    else
        subject.stats = spm_jobman('run',matlabbatch); clear matlabbatch;
    end
end

% task stuff only
% ----------------
if ~exist('resting_state','var')
    
    % contrasts
    % we can assume that 1 contrast per condition is computed, if there
    % is only one run, no need this is the same as betas, if there are
    % seveal runs, simply match condition labels - and do it for each
    % derivatives if any
    if max(size(subject.func)) > 1
        disp('computing contrasts across conditions')
        con_file1 = [Statspath filesep 'con_0001.nii'];
        
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                && ~exist(con_file1,'file'))
            
            SPM        = load([Statspath filesep 'SPM.mat']);
            SPM        = SPM.SPM;
            
            pattern =  'Sn\([0-9]*\) R[0-9]*';
            otherreg   = ~cellfun('isempty', regexp(SPM.xX.name, pattern, 'match'));

            pattern =  'Sn\([0-9]*\) constant*';
            constant   = ~cellfun('isempty', regexp(SPM.xX.name, pattern, 'match'));

            conditions = true(numel(SPM.xX.name), 1);
            conditions(otherreg) = 0;
            conditions(constant) = 0;
            
            SPMnames   = unique(SPM.xX.name(conditions));
            for n=1:length(SPMnames)
                SPMnames{n} = SPMnames{n}(strfind(SPMnames{n},' ')+1:end);
            end
            SPMnames   = unique(SPMnames);
            
            % make contrast vectors
            cindex = 1; c = [];
            for n=1:length(SPMnames)
                test = cell2mat(arrayfun(@(x) contains(x.name,SPMnames{n}),SPM.xX,'UniformOutput',false));
                if sum(test) ~= 0
                    cname{cindex} = SPMnames{n}; %#ok<AGROW>
                    c(cindex,:)   = test; %#ok<AGROW>
                    cindex        = cindex +1;
                end
            end
            
            if ~isempty(c)
                matlabbatch{1}.spm.stats.con.spmmat                             = {SPMmat_file};
                for contrast = size(c,1):-1:1
                    matlabbatch{1}.spm.stats.con.consess{contrast}.tcon.name    = cname{contrast};
                    matlabbatch{1}.spm.stats.con.consess{contrast}.tcon.weights = c(contrast,:);
                    matlabbatch{1}.spm.stats.con.consess{contrast}.tcon.sessrep = 'none';
                end
                matlabbatch{1}.spm.stats.con.delete                             = 1;
                spm_jobman('run',matlabbatch); clear matlabbatch;
            end
        end
    end
    
    %% correct hrf for time dispersion
    if options.derivatives ~= 0 
        outfiles = spmup_hrf_boost(SPMmat_file);
        if isempty(options.skernel)
            options.skernel = [3 3 3];
        end
        subject.boostedhrf = spmup_smooth_boostedfiles(outfiles{1},options.skernel);
        if length(outfiles) == 2
            subject.boostedcontrast = spmup_smooth_boostedfiles(outfiles{2},options.skernel);
        end
    end
    
    %% last QA! carpet plot of residuals
    if strcmp(options.carpet_plot,'on')
        Res = dir([Statspath filesep 'Res_*.nii']);
        for n = size(Res,1):-1:1
            P(n,:) = fullfile(Res(n).folder,Res(n).name);
        end
        % spm_file_merge(spm_vol(P),fullfile(Statspath,'Residuals.nii'))
        
        % create carpet plots
        if ~exist(fullfile(Statspath,'voxplot.fig'), 'file')
            FD = [];
            for frun = 1:size(subject.func)
                FD = [FD; subject.func_qa{frun}.FD]; %#ok<AGROW>
            end
            spmup_timeseriesplot(P, subject.tissues{1},subject.tissues{2},subject.tissues{3}, ...
                'motion',FD,'nuisances','on','correlation','off','clean','off','figure','save');
        end
        
        if strcmp(options.keep_data,'off')
            for n = 1:size(P,1)
                delete(P(n,:));
            end
        end
    end
end



