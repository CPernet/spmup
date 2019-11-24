function outliers = spmup_comp_robust_outliers(time_series)

% Computes robust ouliers of a time series using Carling's k
%
% FORMAT outliers = spmup_comp_robust_outliers(time_series)
%
% INPUT  time_series are the fMRI time courses 
% OUTPUT outliers a binary indicating outliers
%
% The outlier detection is the adjusted Carling's box-plot rule
% ie outside the bound of median+/- k*IQR, 
% k = k=(17.63*n-23.64)/(7.74*n-3.71)
% see Carling, K. (2000). Resistant outlier rules and the non-Gaussian case. 
% Stat. Data Anal. 33, 249–258. 
% <http://www.sciencedirect.com/science/article/pii/S0167947399000572>
% The Inter-Quartile-Range is computed the ideal fourths 
% see D.C. Hoaglin, B. Iglewicz (1987) Fine-tuning some resistant rules for
% outlier labelling. J. Amer. Statist. Assoc., 82 , 1147–1149
% <http://www.tandfonline.com/doi/abs/10.1080/01621459.1986.10478363>
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox


% interquartile range
y = sort(time_series);
j = floor(length(time_series)/4 + 5/12);
g = (length(time_series)/4) - j + (5/12);
k = length(time_series) - j + 1;

ql    = (1-g).*y(j) + g.*y(j+1); % lower quartile
qu    = (1-g).*y(k) + g.*y(k-1); % higher quartile
value = qu-ql;                   % inter-quartile range

% robust outliers
M = median(time_series);
k = (17.63*length(time_series)-23.64)/(7.74*length(time_series)-3.71); % Carling's k
outliers = time_series<(M-k*value) | time_series>(M+k*value);


