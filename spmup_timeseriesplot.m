function M = spmup_timeseriesplot(P,c1,c2,c3,varargin)

% routine to produce plots a la Jonathan Power
%
% FORMAT spmup_timeseriesplot(P,c1, c2, c3, options)
%        M = spmup_timeseriesplot(P,c1, c2, c3, 'motion','on','displacement','on','nuisance','on','design','on')
%
% INPUT P is a cell array for the time series (see spm_select)
%       c1, c2, c3 are the tissue classes derived from the segmentation
%       several options are available
%       'motion' is either 'on' or is a vector associated to the data
%       'displacement' is either on or is a vector associated to the data
%       --> these represent the motion and displacement computed from the
%       motion parameter files
%       'nuisance' is either on or is a matrix (n*2) of noise computed from
%       c2 and c3
%       --> these typically reflect changes in globals, respiration and
%       cardiac cycles that are not well captured by motion correction
%       'design' is the tsv file describing the experiment
%       --> for task fMRI this is useful to check against motion
%       'correlations'  
% OUTPUT a figure (voxplot) showing all gray and white matter voxels
%         in time, associated to the traces in options
%        M is the matrix of voxel by time
%
% Reference: Power, J.D. (2016). A simple but useful way to assess fMRI
% scan qualities. NeuroImage
% <http://dx.doi.org/10.1016/j.neuroimage.2016.08.009>
%       
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

%% check data in
if nargin < 4 
    error('at least 4 input elements expected')
end

V = spm_vol(P);
c1 = spm_read_vols(spm_vol(c1));
c2 = spm_read_vols(spm_vol(c2));
c3 = spm_read_vols(spm_vol(c3));
if any(size(c1)~=V(1).dim) || any(size(c2)~=V(1).dim) || any(size(c3)~=V(1).dim)
    error('segmentation images and data are not of the same dimension')
end

% options
motion = []; displacement = [];
nuisance = []; design = [];
for i=1:length(varargin)
   if strcmpi(varargin{i},'motion') 
       motion = varargin{i+1};
   elseif strcmpi(varargin{i},'displacement') 
       displacement = varargin{i+1};
   elseif strcmpi(varargin{i},'nuisance') 
       nuisance = varargin{i+1};
   elseif strcmpi(varargin{i},'design') 
       design = varargin{i+1};
   end 
end

%% update masks
% make each class mutually exclisive by making sure voxel content is higher
% than the sum of the other 10% 
GM  = c1(c1>(c2+c3)); 
GM  = GM > max(GM)+0.1;
WM  = c2(c2>(c1+c3)); 
WM  = WM > max(WM)+0.1;
CSF = c3(c3>(c1+c2)); 
CSF = CSF > max(CSF)+0.1;
clear c1 c2 c3

%% proceed



