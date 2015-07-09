function CM = connectivity_analysis(varargin)

% This script compute the pair-wise correlations between time series (Images)
% from multiple ROI. It depends on spm to read the data.
%
% FORMAT CM = connectivity_analysis(Images,ROIs,peak_size,peak_image)
%
% INPUT Images is a list of fMRI images names (see spm_select)
%       ROIs is a list of binary masks (see spm_select)
%       lag indicates the possible lag between time series (in TR) -
%           default is 0
%
% OUTPUT CM is the connectivity matrix
%
% examples: CM = connectivity_analysis(myresiduals,myROIs)
%           for each ROI, compute the 1st eigen value across voxels for 
%           each of the 'myresiduals' 
%           CM = connectivity_analysis(myresiduals,myROIs,6,spmT_0001)
%           for each ROI, compute the 1st eigen value across voxels for
%           each of the 'myresiduals' images but from 12mm spheres centred
%           on the max value of the spmT_0001 within the ROIs (allows some
%           variability in the voxel choosen across subjects)
%
% For each ROI, the time series are extracted from each voxel. If peak_size
% is defined, only the X voxels around the max are taken. If peak_image is
% defined this is taken from the peak in that image. Once the time series 
% are obtained, the 1st eigen vector is computed (rather than the mean). 
% Finally, the robust Spearman skipped correlation is computed between each 
% ROI to construct the connectivity matrix CM. You can simply call iteratively 
% this function changing the input Images (and possibly peak_image) to
% obtain 3D matrix of connectivty and assess if the the r values differ
% from 0 or differ between groups. 

%% INPUTS
Images = spm_vol(varargin{1});
ROIs = spm_read_vols(spm_vol(varargin{2}));
x = Images(1).dim(1);
y = Images(1).dim(2);
z = Images(1).dim(3);
nimage = size(Images,1);
[xx,yy,zz,nroi]=size(ROIs);
if x~= xx || y~= yy || z~= zz
    error('images to read data from and ROI images don''t have the same dimensions')
end
clear xx yy zz
lag = 0;

%% LOOP PER ROI

for r=1:nroi
    % message user
    fprintf('extracting data for ROI %g \n',r);
    % find voxels of the ROIs
    [X,Y,Z] = ind2sub([x y z],find(ROIs(:,:,:,r)));
    % get the data
    data = spm_get_data(Images,[X Y Z]');
    % remove NaN if anay
    data(:,isnan(data(1,:))) = [];
    % compute the eigen value
    [m,n]   = size(data);
    if m > n
        [v,s,v] = svd(data'*data);
        s       = diag(s);
        v       = v(:,1);
        u       = data*v/sqrt(s(1));
    else
        [u,s,u] = svd(data*data');
        s       = diag(s);
        u       = u(:,1);
        v       = data'*u/sqrt(s(1));
    end
    d       = sign(sum(v));
    u       = u*d;
    v       = v*d;
    ev(:,r)      = u*sqrt(s(1)/n);
end

%% compute the connectivity matrix

% pairs to test
pairs = nchoosek([1:nroi],2);
CM = eye(nroi);

% test with lag of 1 scan
disp('computing correlations')
for col = 1:size(pairs,1)
    if lag == 0
        C = corr(ev(:,pairs(col,:)));
        CM(pairs(col,1),pairs(col,2)) = C(1,2);
        CM(pairs(col,2),pairs(col,1)) = CM(pairs(col,1),pairs(col,2));
    else
        [C,lag]=xcorr(ev(:,pairs(col,:)),lag,'coeff');
        CM(pairs(col,1),pairs(col,2)) = max(C(:,2));
        CM(pairs(col,2),pairs(col,1)) = CM(pairs(col,1),pairs(col,2));
    end
end
    
