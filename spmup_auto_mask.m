function M = spmup_auto_mask(Y,threshold,fig)

% routine to compute a mask from Y, a 4D time-series of images
% this gives a similar (but more inclusive) mask than SPM. Data are
% smoothed and the average is used as a mask above threshold * by all
% voxels that are non zeros differences in the time series
%
% FORMAT M = spmup_auto_mask(Y,threshold)
%
% INPUT  Y, a 4D time-series of images
%        threshold, the percentage of signal to keep (default = 0.2)
%        fig 'on' or 'off' (default) to image the mask and average image
%
% OUTPUT M the mask image
%
% Cyril Pernet v1 17-Feb-2014
% ----------------------------
% Copyright (c) SPMU+ toolbox

threshold = 0.2;
fig = 'off';
if nargin == 2
    threshold = varargin{2};
elseif nargin == 3
    threshold = varargin{2};
    fig = varagin{3};
end

% compute the mean of normalized smoothed data
% by filling small holes we are more inclusive
for im=1:size(Y,4)
    Ys(im,:,:,:) = smooth3(squeeze(Y(:,:,:,im))); 
end
Ys = (Ys-min(Ys(:)));
Ys = Ys ./ max(Ys(:)); 
Avg = squeeze(mean(Ys,1));

% mask is any non 0 difference between successive scan for average values
% abobve the specified threshold
M = squeeze(any(diff(Ys))) .* squeeze(Avg>threshold);

if strcmp(fig,'on')
    figure('Name','Mask')
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
