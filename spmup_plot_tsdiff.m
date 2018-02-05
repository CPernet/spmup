function h = spmup_plot_tsdiff(timeseries1,timeseries2)

% routine that takes timeseries as input (n*p) matrices and plot the least
% the most different time courses
%
% Cyril Pernet 
% --------------------------------------------------------------------------
% Copyright (C) spmup 2017

A = union(find(sum(timeseries1,1)==0),find(sum(timeseries2,1)==0));
timeseries1(:,A) = []; timeseries2(:,A) = []; % clean up 0s

D = sum((timeseries1 - timeseries2).^2,1); % sum of squared differences
[~,p1] = min(D); % find witch voxel has minimum difference
[~,p2] = max(D); % find witch voxel has maximum difference

h = figure;
set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1); plot(timeseries1(:,p1),'LineWidth',2); hold on; grid on; box on
plot(timeseries2(:,p1),'--','LineWidth',2); title('Least different time series');
subplot(2,1,2); plot(timeseries1(:,p2),'LineWidth',2); hold on; grid on; box on
plot(timeseries2(:,p2),'--','LineWidth',2); title('Most different time series');
