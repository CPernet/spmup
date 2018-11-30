function censoring_regressors = spmup_censoring(data)

% routine that computes robust outliers for each column of the data in
% and return a matrix on cessoring regressors (0s and a 1 per column)
%
% Cyril Pernet
% --------------------------------
% Copyright (C) spmup 2017

[n,p] = size(data);
for column = 1:p
    % interquartile range
    y=sort(data(:,column));
    j=floor(n/4 + 5/12);
    g=(n/4)-j+(5/12);
    ql=(1-g).*y(j)+g.*y(j+1); % lower quartile
    k=n-j+1;
    qu=(1-g).*y(k)+g.*y(k-1); % higher quartile
    value=qu-ql; % inter-quartile range
    
    % robust outliers
    M = nanmedian(data(:,column));
    k=(17.63*n-23.64)/(7.74*n-3.71); % Carling's k
    outliers(:,column)=data(:,column)<(M-k*value) | data(:,column)>(M+k*value);
end

% create the motion censoring regressors simply diagonalizing each outlier
% columns and then removing empty columns
censoring_regressors = [];
for column = 1:p
    censoring_regressors = [censoring_regressors diag(outliers(:,column))];
end
censoring_regressors(:,sum(censoring_regressors,1)==0) = []; % remove empty columns



