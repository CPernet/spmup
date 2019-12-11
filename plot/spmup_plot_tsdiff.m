function h = spmup_plot_tsdiff(varargin)

% routine that takes timeseries images as input , along with coordinates
% and plot the least and most different time courses in the time or
% freqency domain
% 
% FORMAT h = spmup_plot_tsdiff(timeseries1,timeseries2,coordinates,options)
% 
% INPUT: timeseries are the name of images to read (3D or 4D)
%        coordinates are the voxel to read data from [x y z]
%        options are 'time' (default) or 'frequency' 
%
% OUTPUT h is the figure handle
%
% Cyril Pernet 
% --------------------------------------------------------------------------
% Copyright (C) spmup 2018

if nargin<3
    error('At least three inputs expected')
end

% check inputs
assert(exist(varargin{1},'file')==2,'%g doesn''t exist',varargin{1})
assert(exist(varargin{2},'file')==2,'%g doesn''t exist',varargin{2})
if size(varargin{3},2) == 3
    coord = varargin{3}';
else
    coord = varargin{3};
end

% get the data
timeseries1 = spm_get_data(spm_vol(varargin{1}),coord); % get raw time series
timeseries2 = spm_get_data(spm_vol(varargin{2}),coord); 
A = union(find(sum(timeseries1,1)==0),find(sum(timeseries2,1)==0));
timeseries1(:,A) = []; timeseries2(:,A) = []; % clean up 0s

% get power spectrum if frequency as option
if nargin>3 && strcmpi(varargin{4},'frequency')
    timeseries1 = pwelch(timeseries1);
    timeseries2 = pwelch(timeseries2);
end
    
% compute the difference
D = sum((timeseries1 - timeseries2).^2,1); % sum of squared differences
[~,p1] = min(D); % find witch voxel has minimum difference
[~,p2] = max(D); % find witch voxel has maximum difference

% plot the results
h = figure;
try
    img= spm_read_vols(spm_vol(varargin{2}));
catch
    img= spm_read_vols(spm_vol(varargin{1}));
end
img = mean(img,4);

set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
if nargin>3 && strcmpi(varargin{4},'frequency')
    subplot(3,4,[1 2 3 4]); plot(10*log10(timeseries1(:,p1)),'LineWidth',2); hold on; grid on; box on
    plot(10*log10(timeseries2(:,p1)),'--','LineWidth',2); axis tight
    title(sprintf('Least different time series, voxel %g %g %g',coord(1,p1),coord(2,p1),coord(3,p1)));
    subplot(3,4,[9 10 11 12]); plot(10*log10(timeseries1(:,p2)),'LineWidth',2); hold on; grid on; box on
    plot(10*log10(timeseries2(:,p2)),'--','LineWidth',2); axis tight
    title(sprintf('Most different time series, voxel %g %g %g',coord(1,p2),coord(2,p2),coord(3,p2)));
    subplot(3,4,6); data = squeeze(img(:,:,coord(3,p1))); data2 = zeros(size(data));
    data2(coord(1,p2),coord(2,p1)) = max(data(:)); imagesc(data+data2); title('min diff voxel')
    subplot(3,4,7); data = squeeze(img(:,:,coord(3,p2))); data2 = zeros(size(data));
    data2(coord(1,p2),coord(2,p2)) = max(data(:)); imagesc(data+data2); title('max diff voxel')
else
    subplot(3,4,[1 2 3 4]); plot(timeseries1(:,p1),'LineWidth',2); hold on; grid on; box on
    plot(timeseries2(:,p1),'--','LineWidth',2); axis tight
    title(sprintf('Least different time series, voxel %g %g %g',coord(1,p1),coord(2,p1),coord(3,p1)));
    subplot(3,4,[9 10 11 12]); plot(timeseries1(:,p2),'LineWidth',2); hold on; grid on; box on
    plot(timeseries2(:,p2),'--','LineWidth',2);  axis tight
    title(sprintf('Most different time series, voxel %g %g %g',coord(1,p2),coord(2,p2),coord(3,p2)));
    subplot(3,4,6); data = squeeze(img(:,:,coord(3,p1))); data2 = zeros(size(data));
    data2(coord(1,p2),coord(2,p1)) = max(data(:)); imagesc(data+data2); title('min diff voxel')
    subplot(3,4,7); data = squeeze(img(:,:,coord(3,p2))); data2 = zeros(size(data));
    data2(coord(1,p2),coord(2,p2)) = max(data(:)); imagesc(data+data2); title('max diff voxel')
end
