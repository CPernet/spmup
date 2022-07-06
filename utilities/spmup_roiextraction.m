function ROI_data = spmup_roiextraction(varargin)

% routine to extract data for an arbitrary number of ROI and from an
% arbitrary number of images - the routine follows the pattern of spm_voi
%
% FORMAT spmup_roiextraction
%        spmup_roiextraction(ROI,Images,VarName)
%
% INPUT if empty the user is prompted
%       ROI array of ROI (binary) images to read (see spm_select)
%       Images array of fMRI data to read (see spm_select)
%       VarName (optional) is the name of the file to save on drive
%
% OUTPUT the file is saved to the drive - the variable name is ROI_data
%        a structure with n roi fields (ie the number of ROI images in input)
%        with the subfields name. n (the number of voxel in the ROI), and ev
%        ev is the 1st eigen value across voxels, obtained for each Images
%
% Cyril Pernet v1 19-02-2015
% ---------------------------
%  Copyright (C) SPMUP Team 


%% Imputs
if nargin == 0
    ROI    = spm_select(Inf,'image','select ROI images');
    Images = spm_select(Inf,'image','select images to extract data from');
elseif nargin >= 2
    ROI    = varargin{1};
    Images = varargin{2};    
else
    error('no or 2 arguments (ROI and Images) are needed, input error')    
end

%% check size compatibility 
ROI = spm_vol(ROI);
if ~spm_check_orientations(ROI)
    error('some ROI images have different dimensions')
end

Images = spm_vol(Images);
if ~spm_check_orientations(Images)
    error('some of the images to read have different dimensions')
end

%% check ROI and Images match
ROIs = spm_read_vols(ROI); % load each ROI and make a large matrix 
ROIs = ROIs > 0;           % make sure it's binary

[xx,yy,zz,nroi] = size(ROIs);
x               = Images(1).dim(1);
y               = Images(1).dim(2);
z               = Images(1).dim(3);

if x~= xx || y~= yy || z~= zz
    error(sprintf('images to read data from and ROI images don''t have the same dimensions \n you could use spmup_resize to fix ROI images'))
end
clear xx yy zz

%% compute

for r=1:nroi
    [~,ROI_data.roi(r).name] = spm_fileparts(ROI(r).fname);
    fprintf('getting data for %s .. \n',ROI_data.roi(r).name)
    
    % find voxels of the ROIs
    [X,Y,Z] = ind2sub([x y z],find(squeeze(ROIs(:,:,:,r))));
    ROI_data.roi(r).n = length(X);
    
    % get the data (nb images * nb voxel in ROIs)
    data = spm_get_data(Images,[X Y Z]');

    % change NaN to mean
    location = find(isnan(mean(data,2)));
    for l=1:length(location)
        data(location(l),isnan(data(location(l),:))) = nanmean(data(location(l),:));
    end
    
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
    % v     = v*d;
    ROI_data.roi(r).ev = u*sqrt(s(1)/n);
end

%% save
varpath = [];
if nargin == 3
    [varpath,varname] = fileparts(varargin{3});
    if isempty(varpath)
        varpath = fileparts(Images(1).fname); % if only name, save next to the data
    end
end

if nargout == 0 && isempty(varpath) % if no output (i.e. GUI call) and no name input
    varname  = spm_input('save as',1,'s');
    varpath  = uigetdir(pwd,['choose directory to save ' varname '.mat']);
end

save(fullflie(varpath,[varname '.mat']),'ROI_data')

