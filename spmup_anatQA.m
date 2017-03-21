function [anatQA,int_data] = spmup_anatQA(varargin)

% *This implements some of the the anatomical QA measures from the*
% *Preprocessed Connectome Project Quality Assurance Protocol (QAP)*
% <http://preprocessed-connectomes-project.org/quality-assessment-protocol/>
%
% FORMAT: anatQA = spmup_anatQA(anat,c1,c2)
% 
% INPUT: anat is the anatomical image
%        c1 is the gray matter image
%        c2 is the white matter image
%
% OUTPUT anatQA is a structure with the following fields:
%        - SNR : the signal-to-Noise Ratio, ie the mean intensity within gray and white matter divided
%                by the standard deviation of the values outside the brain. Higher values are better.
%        - CNR : the Contrast to Noise Ratio, i.e. the mean of the white matter intensity values minus the mean 
%                of the gray matter intensity values divided by the standard deviation of the values outside the
%                brain. Higher values are better.
%        - FBER: Foreground to Background Energy Ratio, i.e. the variance of voxels in grey and white matter 
%                divided by the variance of voxels outside the brain. Higher values are better.
%        - EFC : Entropy Focus Criterion, i.e. the entropy of voxel intensities proportional to the maximum 
%                possibly entropy for a similarly sized image. Indicates ghosting and head motion-induced blurring. 
%                Lower values are better. See <http://ieeexplore.ieee.org/document/650886/>
%       - asym : a simple measure of assymmetry (Left-Right/Left+Right). Assuming the [0,0,0]
%                coordinate is roughly in the middle of the image (+/- 20 mm), asymmetry is
%                computed for the combined gray/white matter masks - large difference
%                can indicate an issue with the image (for healthy brains)
%
% Note: GM and WM are thresholded by making them mutually exclusive
% The background is found using data outside a large brain mask and
% trimming extreme values
% 
% Cyril Pernet 20 January 2017
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2017

% 2nd input int_data is used for unit testing, this retuns intermediate
% values computed here in a structure 

if nargin ~=3
    error('3 immage names expected as input')
end

AnatV  = spm_vol(varargin{1});
GrayV  = spm_vol(varargin{2});
if any(AnatV.dim ~= GrayV.dim)
    error('Ther anatomical and gray matter images do not have the same size')
end
WhiteV = spm_vol(varargin{3});
if any(AnatV.dim ~= WhiteV.dim)
    error('Ther anatomical and white matter images do not have the same size')
end

%% compute

disp('spmup_anatQA: computing brain mask and background value')

% % define voxels at the edge of the brain 
% x = repmat([1 2 AnatV.dim(1)-1 AnatV.dim(1)],4,1);
% xy = [x(:) repmat([1 2 AnatV.dim(2)-1 AnatV.dim(2)]',4,1)];
% z = repmat([1 2 AnatV.dim(3)-1 AnatV.dim(3)],16,1);
% xyz = [repmat(xy,4,1) z(:)];
% bkgd_threshold = median(spm_get_data(AnatV,xyz'));

% define a large brain mask making sure the background has nothing but background in it
brain_mask = (smooth3(spm_read_vols(GrayV),'box',25)+smooth3(spm_read_vols(WhiteV),'box',25)) > 0.1;
[x,y,z] = ind2sub(AnatV.dim,find(brain_mask==0));
data = sort(spm_get_data(AnatV,[x y z]'));
if isempty(data)
    figure; imagesc(squeeze(brain_mask(:,:,round(AnatV.dim(3)/2)))); title('brain mask')
    error(sprintf('no background data were obtained from the brain mask \n error in image %s\n',AnatV.fname));
end

% compute the noise level using and  50% trimmed data (remove extreme values)
data(isnan(data)) = [];
up = median(data)+iqr(data)/2; down = median(data)-iqr(data)/2;
index = logical((data<up).*(data>down));
std_nonbrain = std(data(index)); 

disp('spmup_anatQA: getting measures')
% make gray/white mutually exclusive
I_GM = spm_read_vols(GrayV) > spm_read_vols(WhiteV);
I_WM = spm_read_vols(WhiteV) > spm_read_vols(GrayV); 

% compute taking voxels in I_GM and I_WM
[x,y,z] = ind2sub(AnatV.dim,find(I_GM));
dataGM = spm_get_data(AnatV,[x y z]');
meanGM = mean(dataGM);

[x,y,z] = ind2sub(AnatV.dim,find(I_WM));
dataWM = spm_get_data(AnatV,[x y z]');
meanWM = mean(dataWM);

anatQA.SNR  = ((meanGM+meanWM)/2)/std_nonbrain;
anatQA.CNR  = (meanWM-meanGM)/std_nonbrain;
anatQA.FBER = var([dataGM dataWM])/std_nonbrain^2;

data = spm_read_vols(AnatV);
data(isnan(data)) = [];
Bmax = sqrt(sum(data(:).^2));
anatQA.EFC = nansum((data(:)./Bmax).*abs(log((data(:)./Bmax))));

% asym
middle_img = abs(AnatV.dim(1).*AnatV.mat(1) / 2);
xoffset = AnatV.mat(1,4);
if xoffset<(middle_img-20) || xoffset>(middle_img+20)
    disp('skipping assymmetry index, the 0 x coordinate doesn''t seem to be in the middle of the image')
else
    voxoffset = abs(xoffset./AnatV.mat(1));
    left = floor(voxoffset-1); right = ceil(voxoffset+1);
    data = (spm_read_vols(GrayV)+spm_read_vols(WhiteV))>0;
    N = data(1:left,:,:)-flipud(data(right:end,:,:)); 
    D = data(1:left,:,:)+flipud(data(right:end,:,:)); 
    asym = N./D; clear N D
    anatQA.asym = mean(asym(data(1:left,:,:)~=0));
end

if nargout == 2
    int_data.meanGM = meanGM;
    int_data.meanWM = meanWM;
    int_data.std_nonbrain = std_nonbrain;
    int_data.varGMWM = var([dataGM dataWM]);
    int_data.Bmax = Bmax;
end

         
         
         
         
         
         

