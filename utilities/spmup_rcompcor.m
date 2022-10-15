function [anoise,tnoise,QC] = spmup_rcompcor(fmridata,brainmask,whitematter,csf,graymatter)

% This function compute the 'compcor' regressors of Behzadi et al. 2007
% <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2214855/>
% One difference from the original paper is that we use more robusts 
% estimates for ROI selection based on median and MAD - the PCA retains
% components above 1% chance under the null
%
% FORMAT: [anoise,tnoise]=rcompcor(fmridata,brainmask,whitematter,csf)
%
% INPUT: fmridata the list of image names for fmri data
%        brainmask the name of the brain mask
%        whitematter the name of the white matter mask 
%        csf the name of the csf mask 
%        graymatter the name of the gray matter mask (only use for QC)
%
% OUTPUT: anoise the compoments related to anatomical ROI
%         tnoise the components related to temporal ROI
%         qc a structure with quality control info and images (can be used
%         by e.g. spm_write_vol if we want to write the images down)
%
% The function relies on SPM <http://www.fil.ion.ucl.ac.uk/spm/> to read
% images and extract data.
%
% Cyril Pernet 10 May 2016
% -------------------------------------------------------------------------

fig = 0; % set to one to visually check what is going on

%% Check data input

if nargin == 0
    [fmridata,sts]= spm_select(Inf,'image','Select fMRI data');
    if sts == 0; disp('selection aborded'); return; end
    [brainmask,sts]= spm_select(1,'image','Select brain mask');
    if sts == 0; disp('selection aborded'); return; end
    [whitematter,sts]= spm_select(1,'image','Select white matter tissue image');
    if sts == 0; disp('selection aborded'); return; end
    [csf,sts]= spm_select(1,'image','Select csf tissue image');
    if sts == 0; disp('selection aborded'); return; end
    if length(nargout) == 3 % if QC requested
        if sts == 0; disp('selection aborded'); return; end
        [graymatter,sts]= spm_select(1,'image','Select gray matter tissue image');
    end
    if sts == 0; disp('selection aborded'); return; end
end

% check fmridata is 4D
Vfmri = spm_vol(fmridata);
if size(Vfmri,1) == 1
    error('fMRI data must be time series')
end
    
% check dimensions
if ischar(brainmask)
    Vbrainmask = spm_vol(brainmask); 
    brainmask = spm_read_vols(Vbrainmask);
end

if ~any(Vfmri(1).dim == size(brainmask))
    error('Dimension issue between mask and fMRI data')
end

if fig == 1
    figure; hold on; title('brain mask')
    for z=1:size(brainmask,3)
        imagesc(brainmask(:,:,z)); pause(0.2)
    end
end

Vgraymatter = spm_vol(graymatter); graymatter = spm_read_vols(Vgraymatter);
if ~any(Vgraymatter.dim == size(brainmask))
    error('Dimension issue between mask and gray matter')
end

Vwhitematter = spm_vol(whitematter); whitematter = spm_read_vols(Vwhitematter);
if ~any(Vwhitematter.dim == size(brainmask))
    error('Dimension issue between mask and white matter')
end

Vcsf = spm_vol(csf); csf = spm_read_vols(Vcsf);
if ~any(Vcsf.dim == size(brainmask))
    error('Dimension issue between mask and CSF')
end
    
    
% check if the masks are binary or not, if not threshold at 0.99 assuming
% these images reflect probability maps (allows removing partial volumes)
graymatter = graymatter.*brainmask;
if length(unique(whitematter)) ~= 2
    graymatter = graymatter > 0.99;
end

whitematter = whitematter.*brainmask;
if length(unique(whitematter)) ~= 2
    whitematter = whitematter > 0.99;
end

csf = csf.*brainmask; length(unique(csf))
if length(unique(whitematter)) ~= 2
    csf = csf > 0.99;
end
    
    
%% anoise
% This noise regressor comes from the white matter and csf images
index      = find(whitematter);
[x,y,z]    = ind2sub(size(whitematter),index);
timeseries = spm_get_data(Vfmri,[x y z]');

index      = find(csf);
[x,y,z]    = ind2sub(size(csf),index);
timeseries = [timeseries spm_get_data(Vfmri,[x y z]')];


%%
% if QC then output the thresholded mask images
if length(nargout) == 3
    [fpath,fname,fext]   = fileparts(Vwhitematter.fname);
    Vwhitematter.fname   = [fpath filesep 'thresholded_' fname fext];
    Vwhitematter.descrip = [Vwhitematter.descrip 'thresholded @ 0.99'];
    QC.whitematter.V     = Vwhitematter;
    QC.whitematter.img   = whitematter;
    
    [fpath,fname,fext]   = fileparts(Vcsf.fname);
    Vcsf.fname           = [fpath filesep 'thresholded_' fname fext];
    Vcsf.descrip         = [Vcsf.descrip 'thresholded @ 0.99'];
    QC.whitematter.V     = Vcsf;
    QC.whitematter.img   = csf;    
end

%%
% get regressors from WM and CSF
disp('getting WM and CSF components')
anoise = decompose(timeseries,fig);

%% tnoise
% Following Lund et al. 2006 <http://www.ncbi.nlm.nih.gov/pubmed/16099175>
% we use areas of high temporal standard deviation to define noise, likely
% reflecting signal from the ventricules, edge regions and vessels. An
% image of the selected voxels is saved allowing post-hoc quality control

%%
% extract time courses from the voxel in the brain mask
index      = find(brainmask);
[x,y,z]    = ind2sub(size(brainmask),index);
timeseries = spm_get_data(Vfmri,[x y z]');
        
%%
% remove low freq from these time courses (linear and quadratic fit)
timeseries = spm_detrend(timeseries,2);

%%
% compute a robust estimate of standard deviation like in Tierney et
% al.2016 <http://www.ncbi.nlm.nih.gov/pubmed/26416652>
% rTSNR = median/median absolute deviate 
SNRimg = NaN(size(brainmask));
rTSNR  = median(timeseries,1) ./ mad(timeseries,1,1);
for v = 1:length(index)
    SNRimg(index(v)) = rTSNR(v);
end

%%
% if QC then output an tSNR image
if length(nargout) == 3
    fpath            = fileparts(Vfmri(1).fname);
    Vfmri(1).fname   = [fpath filesep 'tSRN' fext];
    Vfmri(1).descrip = 'robust temporal SNR image';
    QC.SNR.V         = Vfmri(1);
    QC.SNR.img       = SNRimg;    
end

%%
% sort voxels, keep top 2% in each slice
for zz = min(z):max(z)
    vox            = SNRimg(:,:,zz);
    SNRimg(:,:,zz) = vox.*((vox >= prctile(vox(:),98)));
    if fig == 1
        if zz == min(z)
            figure; hold on; title('high SNR voxels')
        end
        imagesc(brainmask(:,:,zz).*SNRimg(:,:,zz)); pause(0.2)
    end
end
    

%%
% if QC then output an image of selected voxels
if length(nargout) == 3
    QC.SNR.selected_voxels = SNRimg;    
end

%%
% get regressor via SVD
index      = find(SNRimg);
[x,y,z]    = ind2sub(size(SNRimg),index);
timeseries = spm_get_data(Vfmri,[x y z]');
disp('getting high temporal noise components')
tnoise     = decompose(timeseries,0);

end

function compout = decompose(timeseries,fig)

% subfunction to threshold compoments using the broken stick method described
% in Jackson (1993) <http://www.jstor.org/stable/1939574?seq=1#page_scan_tab_contents>
% Instead of analytical solution, we do 1000 Monte-Carlo from normally 
% distributed data from a matrix that have the same rank, giving a 
% distribution of expected principal values - and thus allowing thresholding.

%%
% clean up by removing constant and linear trend column-wise
timeseries = spm_detrend(timeseries,1);

%% 
% decompose:: 1- normalization (same weights per variable) 2-SVD(M*M')
timeseries = zscore(timeseries);  
[V,D]      = eig(timeseries*timeseries');
s          = sort(abs(diag(D)),1,'descend')';
propvar    = s/sum(s);
if fig == 1
    figure; plot(propvar.*100,'o'); hold on
    plot(propvar.*100,'LineWidth',3); xlabel('Eigen Value','FontSize',12);
    ylabel('Proportion of variance explained','FontSize',12); box on; grid on
end

%%
% generate Monte-Carlo data and take the 2.5 top percentille of eigen
% values under the null to threshold the observed data
R  = rank(timeseries);
SS = NaN(R,1000);
spmup_setparallel
parfor MC=1:1000
    SS(:,MC) = sort(abs(eig(randn(R))));
end
PP      = SS./ repmat(sum(SS),[R,1]);
PP      = sort(PP,2); 
H       = sort(PP(:,999),1,'descend'); % keep 1%
keep    = propvar > max(H);
compout = V(:,keep);
if fig == 1
    plot(H.*100,'g','LineWidth',2)
    plot(propvar(keep).*100,'go','LineWidth',3);
end

end