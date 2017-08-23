function outlying_volumes = spmup_outliers(varargin)
%
% this function is similar to 3dToutcount
% http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dToutcount.html
% Note is requires the statistics toolbox (nansum, icdf are called)
%
% FORMAT spmup_despike
%        spmup_despike(P)
%        spmup_despike(P,M)
%
% INPUT if none the user is prompted
%       P the names of the fMRI images (time-series) or the 4D matrix of data
%       M the name of the mask or the 3D binary matrix
%
% OUTPUT 
%
% one looks for outlier voxels based on the Median Absolute Deviation
% Here the MAD is median absolute value of time series minus trend.
% The trend is optained using spn_detrend (2nd order polynomial)
% Outliers are then those voxels with values outside
% k = alphav * sqrt(pi/2) * MAD 
%
% Cyril Pernet 
% --------------------------------------------------------------------------
% Copyright (c) SPM Utility Plus toolbox

if exist('nansum','file') ~= 2
    error('you do not have stats toolbox to perform this operation, sorry')
end

if exist('smooth','file') ~= 2
    error('you need the curve fitting toolbox to perform this operation, sorry')
end

%% check inputs

disp('running spmup_outliers ...')
disp('-------------------------')

%% get data and mask
% memory mapped data
if nargin == 0;
    [P,sts] = spm_select(Inf,'image','select the time series',[],pwd,'.*',1);
    V = spm_vol(P);
    % bypass orientation check allowing to look at raw data
    N = numel(V);
    Y = zeros([V(1).dim(1:3),N]);
    for i=1:N
        for p=1:V(1).dim(3)
            Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
        end
    end
    disp('generating a mask')
    Mask = spmup_auto_mask(V);
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
    else
        if numel(size(P)) == 4 % this is already data in
            Y = P; N = size(Y,4);
        else
            error('input data are not char nor 4D data matrix, please check inputs')
        end
    end
    
    if nargin == 2
        M = varargin{2};
        if ischar(M)
            V = spm_vol(M);
            N = numel(V);
            Mask = zeros([V(1).dim(1:3),N]);
            for i=1:N
                for p=1:V(1).dim(3)
                    Mask(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
                end
            end
        else
            if numel(size(M)) == 4 % this is already data in
                Mask = m; N = size(Mask,4);
            else
                error('input data are not char nor 4D data matrix, please check inputs')
            end
        end
    end
end

disp('looking for outlier volumes')

% detrend, get the MAD and classify
class = NaN(size(Y));
alphav = icdf('Normal',1-(0.001/N),0,1);

for x=1:size(Y,1)
    for y=1:size(Y,2)
        index = find(Mask(x,y,:));
        if ~isempty(index)
            clean_data = spm_detrend(squeeze(Y(x,y,:,:)),2); % detrend
            M = repmat(nanmedian(clean_data,2),1,N); % medians of the time-series
            MAD = nanmedian(abs(clean_data-M),2); % Median absolute deviation of the time series
            k = (alphav * sqrt(pi/2)) .* MAD; % how far is far away
            up = repmat(nanmean(clean_data,2)+k,1,N);
            down = repmat(nanmean(clean_data,2)-k,1,N);
            class(x,y,:,:) = (clean_data > up) + (clean_data < down); % threshold
        end
    end
end

% compute the proportion of outliers per volume
Nb_voxels = nansum(Mask(:));
for im=1:N
    tmp = squeeze(class(:,:,:,im));
    outlying_voxels(im) = (nansum(tmp(:))./Nb_voxels)*100;
end
M = repmat(median(outlying_voxels),1,N);
MAD = median(abs(outlying_voxels-M));
MADN = repmat((MAD./.6745),1,N); % this is almost like 3.5*MAD but better
outlying_volumes = (abs(outlying_voxels-M) ./ MADN) > sqrt(chi2inv(0.975,1));
