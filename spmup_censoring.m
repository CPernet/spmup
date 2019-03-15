function censoring_regressors = spmup_censoring(data)

% routine that computes robust outliers for each column of the data in
% and return a matrix of censsoring regressors (0s and a 1 per column)
%
% FORMAT censoring_regressors = spmup_censoring(data)
%
% INPUT data is a n*m matrix 
% OUTPUT censoring_regressors a matrix with ones in each column for
%        outliers found in coumns of data
%
% Cyril Pernet
% --------------------------------
% Copyright (C) spmup 2017


k = 2.2414; % = sqrt(chi2inv(0.975,1))
[n,p] = size(data);
distance = NaN(n,p);
for p=1:size(data,2)
    tmp = data(:,p);
    points = find(~isnan(tmp));
    tmp(isnan(tmp)) = [];
    
    % compute all distances
    n = length(tmp);
    for i=1:n
        j = points(i);
        indices = [1:n]; indices(i) = [];
        distance(j,p) = median(abs(tmp(i) - tmp(indices)));
    end
    
    % get the S estimator
    % consistency factor c = 1.1926;
    Sn = 1.1926*median(distance(points,p));
    
    % get the outliers in a normal distribution
    outliers(:,p) = (distance(:,p) ./ Sn) > k; % no scaling needed as S estimates already std(data)
    outliers(:,p) = outliers(:,p)+isnan(data(:,p));
end

% create the motion censoring regressors simply diagonalizing each outlier
% columns and then removing empty columns
censoring_regressors = [];
for column = 1:p
    censoring_regressors = [censoring_regressors diag(outliers(:,column))];
end
censoring_regressors(:,sum(censoring_regressors,1)==0) = []; % remove empty columns



