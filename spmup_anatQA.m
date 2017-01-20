function anatQA = spmup_anatQA(varargin)

% *This implements some of the the anatomical QA measures from the*
% *Preprocessed Connectome Project Quality Assurance Protocol (QAP)*
% <http://preprocessed-connectomes-project.org/quality-assessment-protocol/>
%
% FORMAT: spmup_anatQA(anat,c1,c2)
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
% 
% Cyril Pernet 20 January 2017
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2017

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

brain_mask = (smooth3(spm_read_vols(GrayV),'box',25)+smooth3(spm_read_vols(WhiteV),'box',25))>0;
[x,y,z] = ind2sub(AnatV.dim,find(brain_mask==0));
data = sort(spm_get_data(AnatV,[x y z]'));
std_nonbrain = std(data(data<(median(data)+iqr(data)/2))); % in case we have high values (ie brain, ghost) remove top 25%

[x,y,z] = ind2sub(AnatV.dim,find(spm_read_vols(GrayV) > 0.6));
dataGM = spm_get_data(AnatV,[x y z]');
meanGM = mean(data);

[x,y,z] = ind2sub(AnatV.dim,find(spm_read_vols(WhiteV) > 0.6));
dataWM = spm_get_data(AnatV,[x y z]');
meanWM = mean(data);

anatQA.SNR  = ((meanGM+meanWM)/2)/std_nonbrain;
anatQA.CNR  = (meanWM-meanGM)/std_nonbrain;
anatQA.FBER = var([dataGM dataWM])/std_nonbrain^2;

data = spm_read_vols(AnatV);
Bmax = sqrt(sum(data(:).^2));
anatQA.EFC = nansum((data(:)./Bmax).*abs(log((data(:)./Bmax))));



