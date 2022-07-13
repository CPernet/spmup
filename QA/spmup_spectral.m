function [f,outliers] = spmup_spectral(varargin)

% fMRI spectral analysis - slice by slice
% allows checking for artifacts, typically in the coil or in k-space reconstruction
%
% FORMAT [f,outliers] = spmup_spectral
%        [f,outliers] = spmup_spectral(P)
%        [f,outliers] = spmup_spectral(P,M)
%        [f,outliers] = spmup_spectral(P,M,'figure','on/save/off')
%
% INPUT P the fMRI data 
%       M the Mask (optional) to apply
%
% OUTPUT f a matrix of power spectral values [volumes * slice * freq]
%        outliers 0/1 indicates which volumes are outliers, slice by slice
%
% Compute the spatial power spectrum (i.e. abs(fft)^2) per slice using a
% zero-padded 256-by-256 element matrix and then average power values
% of all the same freqencies, returning f. Outliers are computed after 
% averaging all frequencies, computing S-outliers among volumes, for each
% slice.
%
% see also spmup_comp_robust_outliers
%
% Cyril Pernet 
% --------------------------
%  Copyright (C) SPMUP Team 

%% get data and mask
% memory mapped data
if nargin == 0
    [P,sts] = spm_select(Inf,'image' ,'Select your fMRI time series',{},pwd,'.*',Inf);
    if sts == 0
        return
    end
    V = spm_vol(P);
    % bypass orientation check allowing to look at raw data
    N = numel(V);
    Y = zeros([V(1).dim(1:3),N]);
    for i=1:N
        for p=1:V(1).dim(3)
            Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
        end
    end
else
    P = varargin{1};
    if ischar(P)
        V = spm_vol(P);
        N = numel(V);
        Y = zeros([V(1).dim(1:3),N]);
        for i=1:N
            for p=1:V(1).dim(3)
                Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
            end
        end
    elseif iscell(P)
        for v=size(P,1):-1:1
            V(v) =spm_vol(P{v});
        end
        
        for i=numel(V):-1:1
            for p=V(1).dim(3):-1:1
                Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
            end
        end
    else
        if numel(size(P)) == 4 % this is already data in
            Y = P; 
        else
            error('input data are not char nor 4D data matrix, please check inputs')
        end
    end
end
    
if nargin == 2
    M = varargin{2};
    if ischar(M)
        V = spm_vol(M);
        Mask = zeros([V(1).dim(1:3),N]);
        for i=1:numel(V)
            for p=1:V(1).dim(3)
                Mask(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
            end
        end
    else
        if numel(size(M)) == 3 % this is already data in
            Mask = m; 
        else
            error('mask data are not char nor 3D data matrix, please check inputs')
        end
    end
else
    disp('generating a mask')
    Mask = spmup_auto_mask(V);
end

fig = 'save';
if nargin == 3
    fig = 'on';
elseif nargin == 4
    if strcmpi(varargin{3},'figure')
        fig = varargin{4};
    end
end

%% computye the image spectral power slice by slice
nbimage = size(V,1);
f       = zeros(nbimage,V(1).dim(3),256/2+1);
Image   = Y.*repmat(Mask,[1,1,1,nbimage]);
clear   Y

for j=1:nbimage
    for i=1:V(1).dim(3)-1 % for each slice
        im          = squeeze(Image(:,:,i,j)); im(isnan(im)) = 0;
        % figure; subplot(1,3,1); imagesc(im); title('Observed slice')
        N           = min(size(im));
        index       = (max(size(im)) - N) / 2;
        im          = im((1+index):size(im,1)-index,:); % make the image square
        imf2        = fftshift(fft2(im,256,256)); % zero-pads im to be 256-by-256
        freq        = -256/2:256/2-1; 
        impf        = abs(imf2).^2; 
        % subplot(1,3,2); imagesc(freq,freq,log10(impf)); title('slice in Fourier space')
        [X,Y]       = meshgrid(freq,freq); % get coordinates of the power spectrum image
        [~,rho]     = cart2pol(X,Y); % equivalent in polar coordinates
        rho         = round(rho);
        for r =(256/2):-1:0
            f(j,i,r+1) = mean(impf(rho==r)); % average power values of all the same freq
        end
        % subplot(1,3,3); freq2=0:256/2; loglog(freq2,squeeze(f(j,i,:))); title('frequency spectrum','Fontsize',14); axis tight
    end
end

%% figure
avgf     = squeeze(mean(f,1)); % avg over time series
avgv     = squeeze(mean(f,3)); % avg freq to get one value per volume and slice
outliers = spmup_comp_robust_outliers(avgv);

if ~strcmpi(fig,'off') 
    
    figure('Name','Slice spectral analysis');
    if strcmpi(fig,'on')
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    else
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1],'visible','off')
    end
    
    subplot(1,2,1); loglog(0:256/2,avgf);
    title('Slices power spectrum averaged over volumes','LineWidth',2);
    axis tight; grid on;  xlabel('freq'); ylabel('power')
    subplot(1,2,2); semilogy(1:size(avgv,1),avgv); 
    semilogy(1:size(avgv,1),avgv.*outliers,'ro');
    title('Avg power spectra per slice','LineWidth',2);
    axis tight; grid on;  xlabel('volumes'); ylabel('power'); hold on; 
    
    if strcmpi(fig,'save')
        if exist(fullfile(fileparts(V(1).fname),'spmup_QC.ps'),'file')
            print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(fileparts(V(1).fname),'spmup_QC.ps'));
        else
            print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(fileparts(V(1).fname),'spmup_QC.ps'));
        end
        close(gcf)
    end
end

