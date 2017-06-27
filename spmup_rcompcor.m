function [anoise,tnoise]=spmup_rcompcor(fmridata,brainmask,whitematter,csf,graymatter)

% This function compute the 'compcor' regressors of Behzadi et al. 2007
% <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2214855/>
% One difference from the original paper is that we use more robusts 
% estimates for ROI selection
%
% FORMAT: [anoise,tnoise]=spmup_rcompcor(fmridata,brainmask,whitematter,csf,graymatter)
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
% The University of Edinburgh, Neuroimaging sciences
% -------------------------------------------------------------------------

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
end

% check fmridata is 4D
if iscell(fmridata)
    for v=1:size(fmridata,1)
        Vfmri(v) =spm_vol(fmridata{v});
    end
else
    Vfmri = spm_vol(fmridata);
end

if sum(size(Vfmri)) == 2
    error('fMRI data must be time series')
end
    
% check dimensions
Vbrainmask = spm_vol(brainmask); brainmask = spm_read_vols(Vbrainmask);
if ~any(Vfmri(1).dim == Vbrainmask.dim)
    error('Dimension issue between mask and fMRI data')
end

% figure; hold on
% for z=1:Vbrainmask.dim(3)
%     imagesc(brainmask(:,:,z)); pause(0.2)
% end
    
Vwhitematter = spm_vol(whitematter); whitematter = spm_read_vols(Vwhitematter);
if ~any(Vwhitematter.dim == Vbrainmask.dim)
    error('Dimension issue between mask and white matter')
end

Vcsf = spm_vol(csf); csf = spm_read_vols(Vcsf);
if ~any(Vcsf.dim == Vbrainmask.dim)
    error('Dimension issue between mask and CSF')
end

if length(nargout) == 3
    Vgraymatter = spm_vol(graymatter); graymatter = spm_read_vols(Vgraymatter);
    if ~any(Vgraymatter.dim == Vbrainmask.dim)
        error('Dimension issue between mask and gray matter')
    end
end

% check if the masks are binary or not, if not threshold at 0.99 assuming
% these images reflect probability maps (allows removing partial volumes)
whitematter = whitematter.*brainmask;
if length(unique(whitematter)) ~= 2
    whitematter = (whitematter>0.99);
end

csf = csf.*brainmask; 
if length(unique(csf)) ~= 2
    csf = (csf>0.99);
end

if length(nargout) == 3
    graymatter = graymatter.*brainmask;
    if length(unique(whitematter)) ~= 2
        graymatter = graymatter > 0.99;
    end
end
    

%% anoise
% This noise regressor comes from the white matter and csf images
index = find(whitematter);
[x,y,z]=ind2sub(size(whitematter),index);
timeseries = spm_get_data(Vfmri,[x y z]');

index = find(csf);
[x,y,z]=ind2sub(size(csf),index);
timeseries = [timeseries spm_get_data(Vfmri,[x y z]')];

%%
% if QC then output the thresholded mask images
if length(nargout) == 3
    [fpath,fname,fext]=fileparts(Vwhitematter.fname);
    Vwhitematter.fname = [fpath filesep 'thresholded_' fname fext];
    Vwhitematter.descrip = [Vwhitematter.descrip 'thresholded @ 0.99'];
    QC.whitematter.V = Vwhitematter;
    QC.whitematter.img = whitematter;
    
    [fpath,fname,fext]=fileparts(Vcsf.fname);
    Vcsf.fname = [fpath filesep 'thresholded_' fname fext];
    Vcsf.descrip = [Vcsf.descrip 'thresholded @ 0.99'];
    QC.whitematter.V = Vcsf;
    QC.whitematter.img = csf;    
end

%%
% get regressors
anoise = decompose(timeseries,0);


%% tnoise
% Following Lund et al. 2006 <http://www.ncbi.nlm.nih.gov/pubmed/16099175>
% we use areas of high temporal standard deviation to define noise, likely
% reflecting signal from the ventricules, edge regions and vessels. An
% image of the selected voxels is saved allowing post-hoc quality control

%%
% extract time courses from the voxel in the brain mask
index = find(brainmask);
[x,y,z]=ind2sub(size(brainmask),index);
timeseries = spm_get_data(Vfmri,[x y z]');
% figure; plot(mean(timeseries,2))

%%
% remove low freq from these time courses (linear and quadratic fit)
timeseries = spm_detrend(timeseries,2);
% figure; plot(mean(timeseries,2))

%%
% compute a robust estimate of standard deviation like in Tierney et
% al.2016 <http://www.ncbi.nlm.nih.gov/pubmed/26416652>
% rTSNR = median/median absolute deviate 
SNRimg = NaN(size(brainmask));
rTSNR = median(timeseries,1) ./ mad(timeseries,1,1);
for v=1:length(index); SNRimg(index(v)) = rTSNR(v); end
% figure; for z=1:Vbrainmask.dim(3); imagesc(SNRimg(:,:,z)); pause(0.2); end

%%
% if QC then output an tSNR image
if length(nargout) == 3
    [fpath,fname,fext]=fileparts(Vbrainmask.fname);
    Vbrainmask.fname = [fpath filesep 'tSRN' fext];
    Vbrainmask.descrip = ['robust temporal SNR image'];
    QC.SNR.V = Vbrainmask;
    QC.SNR.img = SNRimg;    
end

%%
% sort voxels, keep top 2% in each slice
for z=1:Vbrainmask.dim(3);
    vox = SNRimg(:,:,z); 
    SNRimg(:,:,z) = (vox >= prctile(vox(:),98));
end
% figure; for z=1:Vbrainmask.dim(3); imagesc(SNRimg(:,:,z)); pause(0.2); end

%%
% if QC then output an image of selected voxels
if length(nargout) == 3
    QC.SNR.selected_voxels = SNRimg;    
end

%%
% get regressor via SVD
index = find(SNRimg);
[x,y,z]=ind2sub(size(SNRimg),index);
timeseries = spm_get_data(Vfmri,[x y z]');
% figure; plot(mean(timeseries,2))
tnoise = decompose(timeseries,0);

%% spectral analysis
% normalize the power spectrum of each time voxel time course by the max
% and average for all voxels in the gray matter mask and compare to the
% noise regressors using cross coherence analysis
% if length(nargout) == 3
%     QC.spectral.graymatter = 
%     QC.spectral.whitematter = 
%     QC.spectral.csf= 
%     QC.spectral.gray-white=
%     QC.spectral.gray-csf=
%     QC.spectral.gray-tsnr=
% end

end

function compout = decompose(timeseries,fig)

% subfunction to threshold compoments using the broken stick method described
% in Jackson (1993) <http://www.jstor.org/stable/1939574?seq=1#page_scan_tab_contents>
% Instead of analytical soltion, we do 1000 Monte-Carlo from normally 
% distributed data from a matrix that have the same rank, giving a 
% distribution of expected principal values - and thus allowing thresholding.

%%
% clean up by removing constant and linear trend column-wise
[m,n] = size(timeseries);
timeseries = spm_detrend(timeseries-(ones(m,1)*mean(timeseries)),1);
% figure; plot(mean(timeseries,2))

%% 
% decompose: 1- normalization (same weights per variable) 2-SVD(M*M')
timeseries = timeseries ./ (ones(m,1)*var(timeseries));  
timeseries(:,find(isnan(sum(timeseries,1)))) = [];
% figure; plot(mean(timeseries,2))
[V,D] = eig(timeseries*timeseries');
s = sort(abs(diag(D)),1,'descend')';
propvar = s/sum(s);
if fig ~= 0
    figure; plot(propvar.*100,'o'); hold on
    plot(propvar.*100,'LineWidth',3); xlabel('Eigen Value','FontSize',12);
    ylabel('Proportion of variance explained','FontSize',12); box on; grid on
end

%%
% generate Monte-Carlo data and take the 2.5 top percentille of eigen
% values under the null to threshold the observed data
R = rank(timeseries);
SS = NaN(R,1000);
parfor MC=1:1000
    SS(:,MC) = sort(abs(eig(randn(R))));
end
PP = SS./ repmat(sum(SS),[R,1]);
PP = sort(PP,2); H=sort(PP(:,975),1,'descend');
keep = propvar > max(H);
compout = V(:,keep);
if fig ~= 0
    plot(H.*100,'r','LineWidth',3)
    plot(propvar(keep).*100,'ro','LineWidth',3);
end

end


% V = cov(timeseries); SD = sqrt(diag(V)); [COEFF, LATENT, EXPLAINED] = pcacov(V./(SD*SD'))
