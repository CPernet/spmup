function [outliers] = spmup_comp_robust_outliers(time_series)
% compute robust ouliers of a time series

% interquartile range
y=sort(time_series);
j=floor(length(time_series)/4 + 5/12);
g=(length(time_series)/4) - j + (5/12);
k=length(time_series) - j + 1;

ql=(1-g).*y(j) + g.*y(j+1); % lower quartile
qu=(1-g).*y(k) + g.*y(k-1); % higher quartile
value=qu-ql; % inter-quartile range

% robust outliers
M = median(time_series);
k=(17.63*length(time_series)-23.64)/(7.74*length(time_series)-3.71); % Carling's k
outliers=time_series<(M-k*value) | time_series>(M+k*value);

end

