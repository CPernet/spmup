function [] = spmup_BIDS_1rstlevel(BIDS_dir, BIDS, subjects, s, options)

% level analysis of BIDS fMRI data - with various options available

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% --------------------
% 1st level modelling
% --------------------
filepath = fileparts(fileparts([BIDS_dir filesep subjects{s}.func(frun,:)]));
SPMmat_file = [filepath filesep 'Stats' filesep 'SPM.mat'];
if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') && ~exist(SPMmat_file,'file'))
    
    clear matlabbatch; mkdir([filepath filesep 'Stats'])
    matlabbatch{1}.spm.stats.fmri_spec.dir = {[filepath filesep 'Stats']};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = BIDS.subjects(1).func(1).meta.RepetitionTime;
    N = length(BIDS.subjects(s).func(1).meta.SliceTiming);
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = N;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(N/ 2);
    
    % sessions, onsets and durations
    for frun = 1:size(subjects{s}.func,1)
        [~,filename,~]=fileparts(subjects{s}.func(frun,:));
        [~,~,info_position] = intersect([filename(1:end-4) 'events.tsv'], all_names);
        N_events = length(BIDS.subjects(s).func(info_position).meta.trial_type);
        for n=1:N_events
            cond{n} = cell2mat(BIDS.subjects(s).func(info_position).meta.trial_type(n,:));
        end
        cond = unique(cond);
        all_cond{frun} = cond;
        N_cond = length(cond);
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(frun).scans = {stats_ready{frun}};
        onsets = NaN(N_events,N_cond);
        durations = NaN(N_events,N_cond);
        for C = 1:N_cond
            for n=1:N_events
                if strcmp(cell2mat(BIDS.subjects(s).func(info_position).meta.trial_type(n,:)),cond{C})
                    onsets(n,C) = BIDS.subjects(s).func(info_position).meta.onset(n);
                    durations(n,C) = BIDS.subjects(s).func(info_position).meta.duration(n);
                end
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).name = cond{C};
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).onset = onsets(~isnan(onsets(:,C)),C);
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).duration = durations(~isnan(durations(:,C)),C);
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(frun).cond(C).orth = 1;
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(frun).multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(frun).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(frun).multi_reg = {multi_reg{frun}};
        matlabbatch{1}.spm.stats.fmri_spec.sess(frun).hpf = 128;
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    if strcmp(options.derivatives,'off') % otherwise do it after stats
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    end
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    if BIDS.subjects(1).func(1).meta.RepetitionTime >=1
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    else
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'FAST';
    end
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', ...
        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    save([filepath filesep 'stats_batch_sub' num2str(s) '.mat'],'matlabbatch');
    spm_jobman('run',matlabbatch); clear matlabbatch;
end

% contrasts
% we can assume that at 1 contrast per condition is computed, if there
% is only one run, no need this is the same as betas, if there are
% seveal runs, simply match condition labels - and do it for each
% derivatives if any
if size(subjects{s}.func,1) > 1
    if ~exist('CI','var') %% if we did not recompute the model this is missing
        for frun = 1:size(subjects{s}.func,1)
            [~,filename,~]=fileparts(subjects{s}.func(frun,:));
            [~,~,info_position] = intersect([filename(1:end-4) 'events.tsv'], all_names);
            N_events = length(BIDS.subjects(s).func(info_position).meta.trial_type);
            for n=1:N_events
                cond{n} =cell2mat(BIDS.subjects(s).func(info_position).meta.trial_type(n,:));
            end
            all_cond{frun} = unique(cond);
            if frun > 1
                CI = intersect(all_cond{frun-1},all_cond{frun});
            end
        end
    end
    
    if ~isempty(CI)
        filepath = fileparts(SPMmat_file);
        con_file1 = [filepath filesep 'con_0001.nii'];
        if strcmp(options.overwrite_data,'on') || (strcmp(options.overwrite_data,'off') ...
                && ~exist(con_file1,'file'))
            
            % check contrast weight for run 1 and replicate if needed
            for n=1:length(BIDS.subjects(s).func(info_position).meta.trial_type)
                cond{n} = cell2mat(BIDS.subjects(s).func(info_position).meta.trial_type(n));
            end
            cond = unique(cond);
            v = zeros(length(cond));
            for c=1:length(cond) % that's the columns to span
                for a=1:length(all_cond{1}) % that's the contrast to do
                    if strcmp(all_cond{1}{a},cond{c})
                        v(a,c) = 1;
                    end
                end
            end
            
            matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_file};
            for frun = 1:size(subjects{s}.func,1)
                for c=1:length(CI)
                    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = CI{c};
                    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = v(c,:);
                    if size(BIDS.subjects(s).func,2) == 1
                        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
                    else
                        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'repl';
                    end
                end
            end
            matlabbatch{1}.spm.stats.con.delete = 1;
            spm_jobman('run',matlabbatch); clear matlabbatch;
        end
    end
    
    if strcmp(options.derivatives,'on')
       outfiles = spmup_hrf_boost(SPMmat_file);
       spmup_smooth_boostedfiles(outfiles{1},options.skernel);
       if length(outfiles) == 2
           spmup_smooth_boostedfiles(outfiles{2},options.skernel);
       end
    end
end


end

