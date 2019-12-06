function outliers = spmup_comp_robust_outliers(time_series,type)

% Computes robust ouliers of a time series using S-outliers (default) or Carling's k
%
% FORMAT outliers = spmup_comp_robust_outliers(time_series,type)
%
% INPUTS time_series are the time courses as row vectors 
%        type is 'S-outliers' or 'Carling'
% OUTPUT outliers a binary indicating outliers
%
% S-outliers is the default options, it is independent of a mesure of
% centrality as this is based on the median of pair-wise distances. This is
% a very sensitive measures, i.e. it has a relatively high false positive
% rates. As such it is a great detection tools. 
%
% The adjusted Carling's box-plot rule can also be used, and derived from
% the median of the data: outliers are outside the bound of median+/- k*IQR, 
% k = k=(17.63*n-23.64)/(7.74*n-3.71). This is a more specific measure,
% as such it is 'better' than S-outliers to regress-out, removing bad data
% points (assuming we don;t want to 'remove' too many).
%
% see:
%
% Rousseeuw, P. J., and Croux, C. (1993). Alternatives to the the median 
% absolute deviation. J. Am. Stat. Assoc. 88, 1273–1263.
% <https://www.tandfonline.com/doi/abs/10.1080/01621459.1993.10476408>
% Carling, K. (2000). Resistant outlier rules and the non-Gaussian case. 
% Stat. Data Anal. 33, 249:258. 
% <http://www.sciencedirect.com/science/article/pii/S0167947399000572>
% Hoaglin, D.C., Iglewicz, B. (1987) Fine-tuning some resistant rules for
% outlier labelling. J. Amer. Statist. Assoc., 82 , 1147:1149
% <http://www.tandfonline.com/doi/abs/10.1080/01621459.1986.10478363>
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

if ~exist('type','var')
    type = 'S-outliers';
end

flip = 0; % in case input is a single column vector
if size(time_series,1) == 1
    flip        = 1;
    time_series = time_series';
end

if strcmpi(type,'S-outliers')
    k        = sqrt(chi2inv(0.975,1));
    for p=size(time_series,2):-1:1
        tmp             = time_series(:,p);
        points          = find(~isnan(tmp));
        tmp(isnan(tmp)) = [];
        
        % compute all distances
        n = length(tmp);
        for i=n:-1:1
            j          = points(i);
            indices    = 1:n; 
            indices(i) = []; % all but current data point
            distance(j,p) = median(abs(tmp(i) - tmp(indices))); % median of all pair-wise distances
        end
        
        % get the S estimator
        % consistency factor c = 1.1926;
        Sn = 1.1926*median(distance(points,p));
        
        % get the outliers in a normal distribution
        outliers(:,p) = (distance(:,p) ./ Sn) > k; % no scaling needed as S estimates already std(data)
        outliers(:,p) = outliers(:,p)+isnan(time_series(:,p));
    end
    
elseif strcmpi(type,'Carling')
    
    % interquartile range
    n = size(time_series,1);
    y = sort(time_series);
    j = floor(n/4 + 5/12);
    g = (n/4) - j + (5/12);
    k = n - j + 1;
    
    ql     = (1-g).*y(j,:) + g.*y(j+1,:); % lower quartiles
    qu     = (1-g).*y(k,:) + g.*y(k-1,:); % higher quartiles
    values = qu-ql;                       % inter-quartiles range
    
    % robust outliers
    M = median(time_series);
    k = (17.63*n-23.64)/(7.74*n-3.71); % Carling's k
    outliers = time_series<(M-k*values) | time_series>(M+k*values);

else
    disp('unrecognized outlier type - empty output')
    outliers = [];    
end

if flip == 1
    outliers = outliers';
end
