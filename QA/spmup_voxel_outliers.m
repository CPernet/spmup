function [outlying_volumes,outlying_voxels] = spmup_voxel_outliers(varargin)
%
% this function is similar to 3dToutcount
% http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dToutcount.html
% Note that it requires the statistics toolbox (nansum, icdf are called)
%
% FORMAT [outlying_volumes,outlying_voxels] = spmup_outliers
%        [outlying_volumes,outlying_voxels] = spmup_outliers(P)
%        [outlying_volumes,outlying_voxels] = spmup_outliers(P,M)
%
% INPUT if none the user is prompted
%       P the names of the fMRI images (time-series) or the 4D matrix of data
%       M the mask (optional) to apply
%
% OUTPUT outlying_volumes indicates those volumes which have a high count
%        of outliying voxels (how many voxels per volume)
%
% one looks for outlier voxels based on the Median Absolute Deviation
% Here the MAD is median absolute value of time series minus trend.
% The trend is optained using spm_detrend (2nd order polynomial)
% Outliers are then those voxels with values outside: k = alphav * sqrt(pi/2) * MAD 
% Once obtain we count the number of outlying voxels per volume to find
% volumes with high number of outlying voxels (same procedure)
%
% Cyril Pernet 
% --------------------------------------------------------------------------
% Copyright (c) SPM Utility Plus toolbox

if exist('nansum','file') ~= 2
    error('you do not have stats toolbox to perform this operation, sorry')
end

%% check inputs

disp('running spmup_outliers ...')
disp('-------------------------')

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

%% compute

% detrend, get the MAD and classify
disp('looking for outlier volumes')
class  = NaN(size(Y));
alphav = icdf('Normal',1-(0.001/N),0,1);

for x=1:size(Y,1)
    for y=1:size(Y,2)
        index = find(Mask(x,y,:));
        if ~isempty(index)
            clean_data = spm_detrend(squeeze(Y(x,y,:,:)),2); % detrend
            M          = repmat(nanmedian(clean_data,2),1,N); % medians of the time-series
            MAD        = nanmedian(abs(clean_data-M),2); % Median absolute deviation of the time series
            k          = (alphav * sqrt(pi/2)) .* MAD; % how far is far away
            up         = repmat(nanmean(clean_data,2)+k,1,N);
            down       = repmat(nanmean(clean_data,2)-k,1,N);
            class(x,y,:,:) = (clean_data > up) + (clean_data < down); % threshold
        end
    end
end

% compute the proportion of outliers per volume
Nb_voxels = nansum(Mask(:));
for im=N:-1:N
    tmp                 = squeeze(class(:,:,:,im));
    outlying_voxels(im) = (nansum(tmp(:))./Nb_voxels)*100;
end
outlying_volumes = spmup_comp_robust_outliers(outlying_voxels);
fprintf('outlier volumes: %g\n',find(outlying_volumes));

