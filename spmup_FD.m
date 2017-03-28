function [FD,RMS] = spmup_FD(realignment_file,radius)

% simple routine that computes framewise displacement as a sum (FD) and RMS
% Power et al. (2012) doi:10.1016/j.neuroimage.2011.10.018  
% Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

x = load(realignment_file); % load data
x = x(:,[1:6]); % make sure to pick only the initial 6 if another file is used
x = spm_detrend(x,1); % de-mean and detrend

if ~exist('radius','var')
    radius = 50;
end
x(:,[4 5 6]) = x(:,[4 5 6]).*radius; % angle in radian by average head size = displacement in mm

D=diff(x,1,1); % 1st order derivative
D=[zeros(1,6); D];  % set first row to 0
FD = sum(abs(D),2); % framewise displacement a la Powers
RMS =sqrt(mean(D.^2,2)); % root mean square for each column a la Van Dijk
    