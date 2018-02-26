function M = spmup_auto_mask(varargin)

% Routine to compute a mask from V, a time-series of memory mapped images.
% This gives a similar (but more inclusive) mask than SPM. Data are
% smoothed and the average is used as the mean above threshold * by all
% voxels that have non zeros differences in the time series
%
% FORMAT M = spmup_auto_mask(V,threshold,fig)
%
% INPUT  V memory mapped images (see spm_vol)
%        threshold the percentage of signal to keep (default = 0.2)
%        fig 'on' or 'off' (default) to image the mask and average image
%
% OUTPUT M the mask image
%        if no output, the mask image is writen on the drive next to 
%        the 1st image in V
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

V= varargin{1};
threshold = 0.2;
fig = 'off';
if nargin == 2
    threshold = varargin{2};
elseif nargin == 3
    threshold = varargin{2};
    fig = varargin{3};
end
clear varargin

% compute the mean of normalized smoothed data
% by filling small holes we are more inclusive
if istruct(V)
    if size(V,1) == 1; V=V'; end
    Ys = NaN(size(V,1),V(1).dim(1),V(1).dim(2),V(1).dim(3));
    parfor im=1:size(V,1)
        Ys(im,:,:,:) = smooth3(squeeze(spm_read_vols(V(im))));
    end
elseif iscell(V)
    Ys = NaN(size(V,1),V{1}.dim(1),V{1}.dim(2),V{1}.dim(3));
    parfor im=1:size(V,1)
        Ys(im,:,:,:) = smooth3(squeeze(spm_read_vols(V{im})));
    end
end

Ys = (Ys-nanmin(Ys(:)));
Ys = Ys ./ nanmax(Ys(:)); 
Avg = squeeze(mean(Ys,1));

% mask is any non 0 difference between successive scan for average values
% abobve the specified threshold
M = squeeze(any(diff(Ys))) .* squeeze(Avg>threshold);

if strcmp(fig,'on')
    figure('Name','Mask - press a key to change slice')
    set(gcf,'Color','w')
    colormap('gray')
    for im=1:size(M,3)
        subplot(1,2,1);
        imagesc(flipud(squeeze(M(:,:,im))'));
        title(['slice' num2str(im)],'Fontsize',14);
        axis square; subplot(1,2,2);
        imagesc(flipud(squeeze(Avg(:,:,im))'));
        title('Average smoothed image','Fontsize',14);
        axis square; pause
    end
end

if nargout == 0
       if iscell(V); V = V{1}; end
       [pathstr,~,ext]= fileparts(V(1).fname);
       V(1).fname = [pathstr filesep 'spmup_mask' ext];
       V(1).descrip = 'spmup mask';
       spm_write_vol(V(1),M);
end
