function spmup_plotqc(tablein,fighandle)

% routine to plot data with kernel density
%
% FORMAT spmup_plotqc(tablein,fighandle)
%
% INPUT tablein is a matlab 'table' i.e. a matrix with properties such as labels
%       figurehandle should be either 'new' (default) or 'current' to either
%           create a new figure or plot on whatever current figure is active
%
% example: spmup_plotqc(AnatQA,'new')
%
% Cyril Pernet July 2022
% --------------------------
%  Copyright (C) SPMUP Team 

if spmup_check_matlab_version()
    return
end

%% datain
if nargin == 1
    fighandle = 'new';
else
    if ~strcmpi(fighandle,{'new','current'})
        warning('fighnadle sould be either ''new'' or ''current'', using new as default')
        fighandle = 'new';
    end
end

if ~istable(tablein)
    error('a matlab table in expected as data in')
end

%% plot
if strcmpi(fighandle,'new')
    if ~isempty(tablein.Properties.Description)
        figure('Name',tablein.Properties.Description)
    else
        figure('Name','spmup plot')
    end
    set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
end

% some figure defaults
point_size    = 25;
gp_dispersion = 0.25; 
grouping      = length(tablein.Properties.VariableNames);
gp_index      = linspace(1,1+gp_dispersion*grouping,grouping);
gp_index      = [gp_index max(gp_index)+0.5];
color_scheme  = parula(256);
color_scheme  = color_scheme(round(linspace(1,256,grouping+1)),:);

tiledlayout('flow')
data = tablein.Variables;
for d = 1:grouping
    nexttile
    tmp                         = sort(data(:,d));
    outliers                    = spmup_comp_robust_outliers(tmp,'Carling');
    Y                           = [tmp tmp];  % make it two columns to display
    S                           = repmat(point_size,[length(tmp),1]);  % constant point size
    Y(1:2:length(tmp),1)        = NaN; % plot every other point
    Y(2:2:length(tmp),2)        = NaN;
    outliers                    = single([outliers outliers]); % same for outliers
    outliers(1:2:length(tmp),1) = NaN;
    outliers(2:2:length(tmp),2) = NaN;
    X                           = repmat([0.85 1.15],[length(tmp),1]); % set X values to plot
    for p=1:size(Y,2)
        scatter(X(outliers(:,p)==0,p),Y(outliers(:,p)==0,p),S(outliers(:,p)==0),color_scheme(d,:)); hold on;
        scatter(X(outliers(:,p)==1,p),Y(outliers(:,p)==1,p),S(outliers(:,p)==1).*1.2,[0 0 0]); 
    end
    % add the kernel 
    [bc,K] = ASH(tmp,100);
    bc(K==0)=[]; K(K==0)= [];
    K = (K - min(K)) ./ max(K); % normalize to [0 1] interval
    high=(K/2); low=(-high);
    y1 = plot(high+1.5,bc);% contours
    set(y1,'Color',color_scheme(d,:)); hold on
    if isnumeric(y1)
        y1 = get(y1); % for older matlab versions
    end
    xpoints = [repmat(min(y1.XData),1,length(y1.XData))',y1.XData'];
    filled  = [y1.YData',y1.YData'];
    % check that we have continuous data, otherwise 0 pad
    if diff(xpoints(1,:)) > 0.001*(range(xpoints(:,1)))
        xtmp          = NaN(size(xpoints,1)+1,size(xpoints,2));
        xtmp(1,:)     = [1.5 1.5];
        xtmp(2:end,:) = xpoints;
        xpoints       = xtmp; clear xtmp
        ytmp          = NaN(size(filled,1)+1,size(filled,2));
        ytmp(1,:)     = filled(1,:);
        ytmp(2:end,:) = filled;
        filled        = ytmp; clear ytmp
    end
    
    if diff(xpoints(end,:)) > 0.001*(range(xpoints(:,1)))
        xpoints(end+1,:) = [1.5 1.5];
        filled(end+1,:)  = filled(end,:);
    end
    hold on;
    fillhandle = fill(xpoints,filled,color_scheme(d,:));
    set(fillhandle,'LineWidth',2,'EdgeColor',color_scheme(d,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color

    title(tablein.Properties.VariableNames{d})
    cst = max(abs(diff(tmp))) * 0.1;
    plot_min = min(tmp)-cst;
    plot_max = max(tmp)+cst;
    axis([0.5 2 plot_min plot_max])
end
end

function [bc,K]= ASH(varargin)

% Computes the density estimate of data using ASH 
% - Average Shifted Histogram is taken from Martinez and Martinez - Comp
% statistics handbook with Matlab.
%
% FORMAT [bc,K]=rst_RASH(data,m)
%
% INPUT  data is a vector
%        m is how many histograms to compute (Default = 100);
%
% OUTPUT bc is the bin count
%        K is the estimated histogram / kernel
%

data = varargin{1};
if sum(isnan(data)) ~=0
    warndlg('NaN data included')
end

if size(data,1) == 1
    data = data';
end
n = length(data);
m = 100; % this is the parameter setting the number of hist in ASH and RASH

if nargin == 2
    m = varargin{2};
end
clear varargin

% a typical problem is with data that have a small range like
% probabilities, the meash is not of the right size or contains
% to few bins - the silution here is simply to scale the data
if range(data) <= 1
    data = data.*10;
    scale = 1;
else
    scale = 0;
end

% Average Shifted Histogram
% ----------------------------------
h = 2.15*sqrt(var(data))*n^(-1/5);
delta = h/m;
% 1 make a mesh with size delta
offset = min(diff(data))/2;
if abs(offset) > 1
    offset = 0.5;
end
t0 = min(data) - offset;
tf = max(data) + offset;
nbin = ceil((tf-t0)/delta);
binedge = t0:delta:(t0+delta*nbin);
out = find(binedge>tf);
if out == 1
    binedge(out) = tf;
else
    binedge(out(1)) = tf;
    binedge(out(2:end)) = [];
end
% 2 bin count
nu = histc(data,binedge);
nu = [zeros(1,m-1) nu' zeros(1,m-1)];
% 3 Get the weight vector.
kern = inline('(15/16)*(1-x.^2).^2');
ind = (1-m):(m-1);
den = sum(kern(ind/m));% Get the denominator.
wm = m*(kern(ind/m))/den;% Create the weight vector.
% 4 Get the bin heights over smaller bins.
ASH=zeros(1,nbin);
for k=1:nbin
    ind=k:(2*m+k-2);
    ASH(k)=sum(wm.*nu(ind));
end
K = ASH/(n*h);
bc = t0+((1:nbin)-0.5)*delta;
if scale == 1
    bc = bc./10;
end
end
