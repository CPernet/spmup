function varargout=spmup_rcompcor(varargin)

% This function compute the 'compcor' regressors of Behzadi et al. 2007
% <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2214855/>
% One difference from the original paper is that we use more robusts
% estimates for ROI selection
%
% FORMAT: [anoise,tnoise,QC] = spmup_rcompcor(fmridata,whitematter,csf,graymatter,brainmask);
%         [anoise,tnoise,QC] = spmup_rcompcor(fmridata,whitematter,csf,graymatter); % implicit brain mask
%         [anoise,tnoise]    = spmup_rcompcor(fmridata,whitematter,csf,brainmask); % no QC
%         [anoise,QC]        = spmup_rcompcor(fmridata,whitematter,csf,graymatter); % no temporal noise regressors
%         [tnoise,QC]        = spmup_rcompcor(fmridata,brainmask); % no anatomical noise regressors
%
% INPUT: fmridata              the list of image names for fmri data
%        whitematter           the name of the white matter mask
%        csf                   the name of the csf mask
%        graymatter (optional) the name of the gray matter mask (only use for anat QC)
%        brainmask  (optional) the name of the mask used to derive voxels with hight
%                              temporal variation (if no input it is the sum of previous
%                              entered masks)
%
% OUTPUT: anoise the compoments related to anatomical ROI
%         tnoise the components related to temporal ROI
%         qc     a structure with quality control info and images (can be used
%                by e.g. spm_write_vol if we want to write the images down)
%
% The function relies on SPM <http://www.fil.ion.ucl.ac.uk/spm/> to read
% images and extract data.
%
% Cyril Pernet 10 May 2016; update November 2017
% The University of Edinburgh, Neuroimaging sciences
% -------------------------------------------------------------------------

%% Check data input

if nargin == 0
    [fmridata,sts]= spm_select(Inf,'image','Select fMRI data');
    if sts == 0; disp('selection aborded'); return; end
    [whitematter,sts]= spm_select(1,'image','Select white matter tissue image');
    if sts == 0; disp('selection aborded'); return; end
    [csf,sts]= spm_select(1,'image','Select csf tissue image');
    if sts == 0; disp('selection aborded'); return; end
    [brainmask,sts]= spm_select(1,'image','Select brain mask');
    if sts == 0; disp('selection aborded'); return; end
    if sts == 0; disp('selection aborded'); return; end
    [graymatter,sts]= spm_select(1,'image','Select gray matter tissue image');
end

for inputs = 1:nargin
    % import and check fRI data (all cases)
    if inputs == 1
        if iscell(varargin{1})
            for v=1:size(varargin{1},1)
                Vfmri(v) =spm_vol(varargin{1}{v});
            end
        else
            Vfmri = spm_vol(varargin{1});
        end
    end
    
    if sum(size(Vfmri)) == 2
        error('fMRI data must be time series')
    end
    
    % case we do it all
    if nargout == 3 && nargin == 5
        if inputs == 1; anoise = []; tnoise = []; QC = []; end
        if inputs == 2; [Vwhitematter,whitematter]=local_import(varargin{2}); end
        if inputs == 3; [Vcsf        ,csf        ]=local_import(varargin{3}); end
        if inputs == 4; [Vgraymatter ,graymatter ]=local_import(varargin{4}); end
        if inputs == 5; [Vbrainmask  ,brainmask  ]=local_import(varargin{5}); end
    end
    
    % case we have no QC
    if nargout == 3 && nargin == 4
        if inputs == 1; anoise = []; tnoise = []; end
        if inputs == 2; [Vwhitematter,whitematter]=local_import(varargin{2}); end
        if inputs == 3; [Vcsf        ,csf        ]=local_import(varargin{3}); end
        if inputs == 4; [Vbrainmask  ,brainmask  ]=local_import(varargin{4}); end
    end
    
    % case we have anatomical noise regressor only
    if nargout <= 2 && nargin == 4
        if inputs == 1 && nargout == 1
            anoise = [];
        else
            anoise = []; QC = [];
        end
        if inputs == 2; [Vwhitematter,whitematter]=local_import(varargin{2}); end
        if inputs == 3; [Vcsf        ,csf        ]=local_import(varargin{3}); end
        if inputs == 4; [Vgraymatter ,graymatter ]=local_import(varargin{4}); end
    end
    
    % case we have temporal noise regressor only
    if nargout <= 2 && nargin == 2
        if inputs == 1 && nargout == 1
            tnoise = [];
        else
            tnoise = []; QC = [];
        end
        if inputs == 2; [Vbrainmask  ,brainmask  ]=local_import(varargin{2}); end
    end
end

if nargout == 3
    QC = [];
end

% quickly check dimensions and if the masks are binary.
% If masks are not binary, threshold at 5% of the max assuming these images
% reflect probability maps (allows removing partial volumes)

if exist('whitematter','var')
    if ~any(Vwhitematter.dim == Vfmri(1).dim)
        error('Dimension issue between mask and white matter')
    end
end

if exist('csf','var')
    if ~any(Vcsf.dim == Vfmri(1).dim)
        error('Dimension issue between mask and gray matter')
    end
end

if exist('graymatter','var')
    if ~any(Vgraymatter.dim == Vfmri(1).dim)
        error('Dimension issue between mask and fMRI data')
    end
end

% we always need a brain mask, if not present, create it
if exist('brainmask','var')
    if ~any(Vbrainmask.dim == Vfmri(1).dim)
        error('Dimension issue between mask and CSF')
    end
else
    if exist('graymatter','var')
        brainmask = whitematter+csf+graymatter;
    else
        brainmask = whitematter+csf;
    end
end
brainmask = brainmask>0;

% figure; hold on
% for z=1:size(brainmask,3)
%     imagesc(brainmask(:,:,z)); pause(0.2)
% end

% make sure other masks are binary
if exist('whitematter','var')
    whitematter = whitematter.*brainmask;
    if length(unique(whitematter)) ~= 2
        values = unique(whitematter);
        whitematter = (whitematter>=values(round(length(values)*0.95)));
    end
end

if exist('csf','var')
    csf = csf.*brainmask;
    if length(unique(csf)) ~= 2
        values = unique(csf);
        csf = (csf>=values(round(length(values)*0.95)));
    end
end

if exist('graymatter','var')
    graymatter = graymatter.*brainmask;
    if length(unique(whitematter)) ~= 2
        values = unique(graymatter);
        graymatter = (graymatter >=values(round(length(values)*0.95)));
    end
end

% ------------------------------------------------------------
%%                   anatomical noise
% ------------------------------------------------------------

if exist('anoise','var')
    % This noise regressor comes from the white matter and csf images
    index = find(whitematter);
    [x,y,z]=ind2sub(size(whitematter),index);
    timeseries = spm_get_data(Vfmri,[x y z]');
    
    index = find(csf);
    [x,y,z]=ind2sub(size(csf),index);
    timeseries = [timeseries spm_get_data(Vfmri,[x y z]')];
    
    % get regressors
    anoise = decompose(timeseries,0);
    
    % if QC then output the thresholded mask images
    if exist('QC','var')
        [fpath,fname,fext]=fileparts(Vwhitematter.fname);
        Vwhitematter.fname = [fpath filesep 'thresholded_' fname fext];
        Vwhitematter.descrip = [Vwhitematter.descrip 'thresholded @ 5% of the max'];
        QC.whitematter.V = Vwhitematter;
        QC.whitematter.img = whitematter;
        
        [fpath,fname,fext]=fileparts(Vcsf.fname);
        Vcsf.fname = [fpath filesep 'thresholded_' fname fext];
        Vcsf.descrip = [Vcsf.descrip 'thresholded @ 5% of the max'];
        QC.csf.V = Vcsf;
        QC.csf.img = csf;
        
        %% spectral analysis
        % normalize the power spectrum of each time voxel time course by the max
        % and average for all voxels in the gray matter mask and compare to the
        % noise regressors using cross coherence analysis
        
        whitematter = timeseries(:,1:end-length(x));
        csf = timeseries(:,end-length(x)+1:end);
        QC.spectral.whitematter = norm_power(whitematter);
        QC.spectral.csf         = norm_power(csf);
        QC.spectral.anoise      = norm_power(anoise);
        
        if exist('graymatter','var')
            [x,y,z]=ind2sub(size(graymatter),find(graymatter));
            QC.spectral.graymatter  = norm_power(spm_get_data(Vfmri,[x y z]'));
        end
        
        figure('Name','Power Spectrum')
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
        plot(QC.spectral.whitematter,'LineWidth',2); hold on
        plot(QC.spectral.csf,'LineWidth',2)
        plot(QC.spectral.anoise,'LineWidth',2)
        if exist('graymatter','var')
            plot(QC.spectral.graymatter,'LineWidth',2)
            legend('WM','CSF','Anat noise','GM','Location','northeast')
            title(['Correlation GM/anat noise:' num2str(max(xcorr(QC.spectral.graymatter,QC.spectral.anoise,'coeff')))])
        else
            legend('WM','CSF','Anat noise','Location','northeast')
            title('WM CSF and Anat noise power spectra')
        end
        grid on; box on; ylabel('power'); xlabel('frequency')
        try
            print (gcf,'-dpsc2', '-bestfit', [pwd filesep 'Anatomical_noise_PCA_regressor.ps']);
        catch
            print (gcf,'-dpsc2', [pwd filesep 'Anatomical_noise_PCA_regressor.ps']);
        end
        close(gcf)
    end
    clear Vwhitematter whitematter Vcsf csf timeseries
    if exist('graymatter','var')
        clear Vgraymatter graymatter
    end
    
    
    varargout{1} = anoise;
    if nargout == 2
        varargout{2} = QC;
    elseif nargout == 3
        varargout{3} = QC;
    end
end

% ------------------------------------------------------------
%%                        temporal noise
% ------------------------------------------------------------
% Following Lund et al. 2006 <http://www.ncbi.nlm.nih.gov/pubmed/16099175>
% we use areas of high temporal standard deviation to define noise, likely
% reflecting signal from the ventricules, edge regions and vessels. An
% image of the selected voxels is saved allowing post-hoc quality control

% extract time courses from the voxel in the brain mask
if exist('tnoise','var')
    
    index = find(brainmask);
    [x,y,z]=ind2sub(size(brainmask),index);
    timeseries = spm_get_data(Vfmri,[x y z]');
    % figure; plot(mean(timeseries,2))
    
    % remove low freq from these time courses (linear and quadratic fit)
    timeseries = spm_detrend(timeseries,2);
    % figure; plot(mean(timeseries,2))
    
    % compute a robust estimate of standard deviation like in Tierney et
    % al.2016 <http://www.ncbi.nlm.nih.gov/pubmed/26416652>
    % rTSNR = median/median absolute deviate
    SNRimg = NaN(size(brainmask));
    rTSNR = median(timeseries,1) ./ mad(timeseries,1,1);
    for v=1:length(index)
        SNRimg(index(v)) = rTSNR(v);
    end
    % figure; for z=1:Vbrainmask.dim(3); imagesc(SNRimg(:,:,z)); pause(0.2); end
    
    % sort voxels, keep top 2% in each slice
    for z=1:Vbrainmask.dim(3)
        vox = SNRimg(:,:,z);
        SNRimg(:,:,z) = (vox >= prctile(vox(:),98));
    end
    % figure; for z=1:Vbrainmask.dim(3); imagesc(SNRimg(:,:,z)); pause(0.2); end
    
    % get regressor via SVD
    index = find(SNRimg);
    [x,y,z]=ind2sub(size(SNRimg),index);
    timeseries = spm_get_data(Vfmri,[x y z]');
    % figure; plot(mean(timeseries,2))
    tnoise = decompose(timeseries,0);
    
    % if QC then output an tSNR image and an image of selected voxels
    if exist('QC','var')
        [fpath,~,fext]=fileparts(Vbrainmask.fname);
        Vbrainmask.fname = [fpath filesep 'tSRN' fext];
        Vbrainmask.descrip = 'robust temporal SNR image';
        QC.SNR.V = Vbrainmask;
        QC.SNR.img = SNRimg;
        QC.spectral.tnoise = norm_power(tnoise);
    end
    
    if nargout == 3 && nargin >= 4
        varargout{2} = tnoise;
        varargout{3} = QC;
    elseif nargout <= 2
        varargout{1} = tnoise;
        if exist('QC','var')
            varargout{2} = QC;
        end
    end
    
end

end

function [V,data]=local_import(name)
V = spm_vol(name);
data = spm_read_vols(V);
end

function compout = decompose(timeseries,fig)

% subfunction to threshold compoments using the broken stick method described
% in Jackson (1993) <http://www.jstor.org/stable/1939574?seq=1#page_scan_tab_contents>
% Instead of analytical soltion, we do 1000 Monte-Carlo from normally
% distributed data from a matrix that have the same rank, giving a
% distribution of expected principal values - and thus allowing thresholding.

%%
% clean up by removing constant and linear trend column-wise
[m,~] = size(timeseries);
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

% V = cov(timeseries); SD = sqrt(diag(V)); [COEFF, LATENT, EXPLAINED] = pcacov(V./(SD*SD'))
end


function AvgP = norm_power(Data)
% routine to return the average power spectrum among voxels (normalized by the max)

[p,~] = size(Data); % how many data point
nfft  = 2^nextpow2(2*p-1);
Data  = spm_detrend(Data,1);
P     = fft(Data,nfft).*conj(fft(Data,nfft)); % compute power
P = P(1:p,:)'; % take only half of the window
AvgP  = nanmean(P./repmat(max(P,[],2),[1 size(P,2)]),1); % normalize each voxel by it's max value in time (P is positive) and average
end

