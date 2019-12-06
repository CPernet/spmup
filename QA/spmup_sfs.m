function sfs = spmup_sfs(fmridata,TR,csf,roi)

% computes the Signal Fluctuation Sensitivity metric
% DeDora et 2016 <http://journal.frontiersin.org/article/10.3389/fnins.2016.00180/full>
%
% FORMAT: SFS = spmup_sfs(fmridata,TR,csf)
%         SFS = spmup_sfs(fmridata,TR,csf,roi)
%
% INPUT: fmridata is the name of the time series (3D or 4D) with has been
%        realigned (and normalized if one wants to use the default roi)
%        scan_param is a structure with fields .TR and .nbslice and coresponding values 
%        csf is the mask for CSF (if probabilistic thresholded at 70%)
%        roi is the binary region of interest to analyze (coud be a network, the
%        computation are simply done for all on zeros voxels)
%        if roi is not in input, the 7 Network of Yeo et al 2011 are analyzed
% 
% OUTPUT: SFS is the Signal Fluctuation Sensitivity defined as
%         mu_roi/mu_global * std_roi/std_csf
%         mu_roi and std_roi are the mean signal and mean std of voxels defined by the ROI image
%         mu_global: the mean signal in time over all in brain voxels
%         std_csf: the mean std in time of CSF voxels (defined as mask > 70%)
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

current = pwd;

%% inputs
if nargin < 3
    fmridata = spm_select(Inf,'image','select fMRI data',[],pwd);
    csf      = spm_select(Inf,'image','select csf mask',[],pwd);
    if isempty(fmridata) || isempty(csf)
        return
    end
    
    if size(fmridata,1) == 1 && strcmp(fmridata(length(fmridata)-1:end),',1')
        fmridata = fmridata(1:length(fmridata)-2); % in case picked 4D put left ,1
    end
    
    TR = input('input TR: ');
    if isempty(TR)
        return
    end
    roi = [fileparts(fileparts(which('spmup_sfs.m'))) filesep 'external' filesep 'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'];
elseif nargin < 4
    % if no ROI specifiy pick up Yeo et al 2011 networks
    roi = [fileparts(fileparts(which('spmup_sfs.m'))) filesep 'external' filesep 'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'];
end

Vfmri = spm_vol(fmridata);
Vcsf = spm_vol(csf); 
if any(Vcsf.dim ~= Vfmri(1).dim)
    error('Dimension issue between data and CSF')
end

% check if the masks are binary with a threshold at 0.7
csf = spm_read_vols(Vcsf);
csf(find(csf<0.7))==0;
csf(find(csf))=1;

%% 1st data scrubbing without filtering
spm('Defaults','fmri')
spm_jobman('initcfg');
out = spmup_realign_qa(fmridata,'Motion Parameters','on',...
    'Voltera','on','Framewise displacement','on','Globals','off','Movie','off');
out = out{2}; % ie the augmented design matrix
mkdir([fileparts(out) filesep 'tmp_sfs'])
matlabbatch{1}.spm.stats.fmri_spec.dir = {[fileparts(out) filesep 'tmp_sfs']};
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {out};
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cell(size(Vfmri,1),1);
for f=1:size(Vfmri,1)
    if size(fmridata,1) == 1
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans{f} = [fmridata ',' num2str(f)];
    else
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans{f} = deblank(fmridata(f,:));
    end
end
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t =  Vfmri(1).dim(3);
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round( Vfmri(1).dim(3)/2);
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = Inf; % infinite high pass
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'none'; % noi low pass
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
newdata = spm_jobman('run',matlabbatch);
Vfmri = spm_vol(cell2mat(newdata{2}.res));

%% 2nd get brain mask and average over voxels the mean signal time course
% from the residuals (data detrended, motion corrected and scrubbed)
disp('estimating brain mask from denoised images and mean signal')
mask = spmup_auto_mask(Vfmri);
[x,y,z]=ind2sub(size(mask),find(mask));
tmp = spm_get_data(Vfmri,[x y z]');
tmp(:,find(sum(isnan(tmp)))) = []; % removes columns of NaN;
mu_global = mean(mean(tmp,1)); % average in time then over voxels

%% 3rd get std from CSF
[x,y,z]=ind2sub(size(csf),find(csf));
tmp = spm_get_data(Vfmri,[x y z]');
tmp(:,find(sum(isnan(tmp)))) = []; % removes columns of NaN;
std_csf = mean(std(tmp,0,1));

%% 4th get ROI values
Vroi = spm_vol(roi);
if any(Vroi.dim ~= Vfmri(1).dim) % maybe we need to resize the ROI image
    clear matlabbatch
    matlabbatch{1}.spm.util.bbox.image = {newdata{2}.res{1}};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    bb = spm_jobman('run',matlabbatch);
    Vroi = spm_vol(cell2mat(spmup_resize(roi,bb{1}.bb,abs(diag(Vfmri(1).mat(1:3,1:3)))')));
end

ROIdata = round(spm_read_vols(Vroi));
roi_values = unique(ROIdata);
roi_values(roi_values==0) = [];
for r = 1:length(roi_values)
    [x,y,z]=ind2sub(size(ROIdata),find(ROIdata==r));
    tmp = spm_get_data(Vfmri,[x y z]');
    tmp(:,find(sum(isnan(tmp)))) = []; % removes columns of NaN;
    sfs(r) = mean(mean(tmp,1))/mu_global * mean(std(tmp,0,1))/std_csf;
end

cd(current)
rmdir('tmp_sfs','s')
