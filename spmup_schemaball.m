function spmup_schemaball(Data,Mask,varargin)

% Routione to plot connectivity matrices
%
% FORMAT spmup_schemaball(Data,Mask,Names,'Threshold',x,'Reorder','Yes',ExtraPlots,'Yes'))
%
% IMPUT Data is a 2D n*n symmetric correlation matrix
%       Mask is a 2D n*n symetric binary correlation matrix
%       OPTIONS ARE
%       Names cell array with n names (default 1 to n)
%       Edges 'Bezier' (default) or 'Line' for curved or straight lines between nodes
%       Connect 'Thickness' or 'Colour' edges change thickness or color based on the Data values
%       Reorder reorder data to minimize edge crossing (bandwidth)
%       ExtraPlots in addition to the circular plot, imagesc of the full and sparsified matrices
%
% Example
% Names{1} = 'Left ROI1'
% Names{2} = 'Left ROI2'
% Names{3} = 'Left ROI3'
% Names{4} = 'Left ROI4'
% Names{5} = 'Right ROI1'
% Names{6} = 'Right ROI2'
% Names{7} = 'Right ROI3'
% Names{8} = 'Right ROI4'
% Data = zeros(8,8); Data(find(eye(8))) = 1;
% Data(1,5) = 0.9; Data(2,6) = 0.8; Data(3,7) = 0.9; % interhemispheric +connections
% Data(1,2) = 0.6; Data(1,3) = 0.6; Data(1,4) = 0.5 % connection ROI1 to others
% Data(5,6) = -0.6; Data(5,7) = -0.6; Data(5,8) = -0.5
% Data = triu(Data).'+triu(Data,1); % = symmetric
% Mask = Data; schemaball(Data,Mask,'Names',Names,'Edges','Bezier','Connect','Colour','Reorder','No','ExtraPlots','Yes')
%
% This version of schemaball is a hack of functions from
% Gunther Struyf https://github.com/GuntherStruyf/matlab-tools/blob/master/schemaball.m
% Oleg Komarov http://www.mathworks.co.uk/matlabcentral/fileexchange/42279-schemaball
%
% Cyril Pernet 13 June 2014

%% check inputs

% set defaults
for i=1:size(Data,1); Names{i}= i; end;
flag.Edges = 'Bezier'; % default is to plot usisng curves
flag.Connect = 'Colour'; % default is to change colour with connection strengh
flag.Reorder = 'Yes'; % default is to reorder data avoiding 0 near diagonal
flag.ExtraPlots = 'Yes'; % default is to also plot the sparsity before/after ordering and the correlation matrix
narginchk(2,12)

for arg = 1:2:(nargin-2)
    if strcmp(varargin(arg),'Names')
        clear Names;
        for n=1:size(varargin{arg+1},2)
            Names{n} = varargin{arg+1}{n};
        end
    end
    
    if strcmp(varargin(arg),'Egdes')
        if strcmp(varargin(arg+1),'Bezier') || strcmp(varargin(arg+1),'Line')
            flag.Edges = varargin(arg+1);
        else
            error('Edge argument must be Bezier or Line')
        end
    end
    
    if strcmp(varargin(arg),'Connect')
        if strcmp(varargin(arg+1),'Thickness') || strcmp(varargin(arg+1),'Colour')
            flag.Connect = varargin(arg+1);
        else
            error('Connect argument must be Thickness or Colour')
        end
    end
    
    if strcmp(varargin(arg),'Reorder')
        if strcmp(varargin(arg+1),'Yes') || strcmp(varargin(arg+1),'No')
            flag.Reorder = varargin(arg+1);
        else
            error('Reorder argument must be Yes or No')
        end
    end
    
    if strcmp(varargin(arg),'ExtraPlots')
        if strcmp(varargin(arg+1),'Yes') || strcmp(varargin(arg+1),'No')
            flag.ExtraPlots = varargin(arg+1);
        else
            error('ExtraPlots argument must be Yes or No')
        end
    end
end

[a,b]=size(Data);
if a~=b
    error('the correlation matrix need top be squared')
else
    N = a; clear a b
end

if size(Data) ~= size(Mask)
    error('The correrlation and binary matrices must be of same size')
end

if length(Names) ~= N
    error('the number of nodes and the number of names differ')
end

% look at the data range to set thickness (from 2 to 6) or colour
if strcmp(flag.Connect,'Thickness')
    n = length(find(tril(Mask,-1)));
    T = exp(1):exp(2.5)/(n+2):exp(2.5);
else
    C = cubehelixmap('continuous',N);
end

%% prepare information
sparsity = nnz(Mask)/numel(Mask);
% path_length = Mask^2;
if strcmp(flag.Reorder,'Yes')
    p = symrcm(Mask); % find permutations to put zeros close to diagonal
    Data2 = Data(p,p);
    Mask2 = Mask(p,p);
    Names2 = Names(p);
    if strcmp(flag.ExtraPlots,'Yes')
        figure; set(gcf,'Color','w');
        subplot(1,2,1); spy(Mask,30); axis([0.5 N+0.5 0.5 N+0.5]);
        title(['Data - sparsity ' num2str(sparsity)])
        set(gca,'XTickLabel',Names); set(gca,'YTickLabel',Names);
        subplot(1,2,2); spy(Mask2,30); axis([0.5 N+0.5 0.5 N+0.5]);
        title('Data with reduced bandwidth');
        set(gca,'XTickLabel',Names2); set(gca,'YTickLabel',Names2);
    end
else
    Data2 = Data;
    Mask2 = Mask;
    Names2 = Names;
end

%% plot a simple correlation matrix
% ------------------------------
if strcmp(flag.ExtraPlots,'Yes')
    figure; set(gcf,'Color','w'); 
    imagesc(mean(Data,3)); colormap(C);colorbar('EastOutside');
    set(gca,'XTickLabel',Names); set(gca,'YTickLabel',Names)
    title('Correlation Matrix','FontSize',16)
end

%% do the schemaball
% -----------------
r = 1;
theta = linspace(0,2*pi,N+1)'; theta(end) = [];
xy = r .* [cos(theta) sin(theta)]; % gives N coordinates
Mask2 = tril(Mask2);
clear Mask Data

% show nodes
figure; set(gcf,'Color','w');
line(xy(:,1), xy(:,2), 'LineStyle','none', 'Marker','o', ...
    'MarkerSize',11, 'MarkerFaceColor',[1 0 0], 'Color','r')
title('Connectivity graph','FontSize',16); hold on
h = text(xy(:,1).*1.05, xy(:,2).*1.05, Names2, 'FontSize',8);

% orientation = -theta*180/pi;
% orientation(orientation>-180) = abs(orientation(orientation>-180));
% set(h, {'Rotation'},num2cell(orientation))
% show edges
if strcmp(flag.Edges, 'Line')
    % ideally we would use color gradient based on Data2
    % figure properties are not easy to find out - help?
    % set lineWidth to 3
    gplot(Mask2.*Data2>0, xy, 'b');
    gplot(Mask2.*Data2<0, xy, 'r');
    % if strcmp(flag.Connect,'Thickness')
    %
    % else
    %
    % end
    axis equal off
elseif strcmp(flag.Edges, 'Bezier');
    % compute Bezier curves
    % Calculate Bx and By positions of quadratic Bezier curves with P1 at (0,0)
    % B(t) = (1-t)^2*P0 + t^2*P2 where t is a vector of points in [0, 1] and determines, i.e.
    % how many points are used for each curve, and P0-P2 is the node pair with (x,y) coordinates.
    % loop per Edge for negative correlations
    tmp = Mask2.*(Data2<0);
    [ind1,ind2]=ind2sub(size(Mask2),find(tmp(:))); % look for non-zeros
    N2 = length(ind1); % number of edges to plot
    colortick1 = [0:N2];
    t = (0.025: 0.05 :1)';
    t2 = [1-t, t].^2;
    Bx = NaN(length(t2),N2);
    By = NaN(length(t2),N2);
    for c = 1:N2
        Bx(:,c) = t2*[xy(ind1(c),1); xy(ind2(c),1)];
        By(:,c) = t2*[xy(ind1(c),2); xy(ind2(c),2)];
    end
    
    % plot but check the correlation strengh in Data2 to choose color order
    cm = cubehelixmap('decrease',N2);
    
    for c = 1:N2
        if strcmp(flag.Connect,'Thickness')
            plot(Bx(:,c),By(:,c),'b','LineWidth',T(plotindex));
        else
            plot(Bx(:,c),By(:,c),'Color',cm(c,:),'LineWidth',2);
        end
    end
    clear tmp ind1 ind2 N2 Bx By
    
    % -------------------------------------
    % loop per Edge for positive correlations
    tmp = Mask2.*(Data2>0);
    [ind1,ind2]=ind2sub(size(Mask2),find(tmp(:))); % look for non-zeros
    N2 = length(ind1); % number of edges to plot
    colortick2 = [colortick1+1:N2];
    t = (0.025: 0.05 :1)';
    t2 = [1-t, t].^2;
    Bx = NaN(length(t2),N2);
    By = NaN(length(t2),N2);
    for c = 1:N2
        Bx(:,c) = t2*[xy(ind1(c),1); xy(ind2(c),1)];
        By(:,c) = t2*[xy(ind1(c),2); xy(ind2(c),2)];
    end
    
    % plot but check the correlation strengh in Data2 to choose color order
     cm = cubehelixmap('increase',N2);
    for c = 1:N2
        if strcmp(flag.Connect,'Thickness')
            plot(Bx(:,c),By(:,c),'b','LineWidth',T(plotindex));
        else
            plot(Bx(:,c),By(:,c),'Color',cm(c,:),'LineWidth',2);
        end
    end
    clear tmp ind1 ind2 N2 Bx By
else
    error('unrecognized Ege option')
end
% Axis settings
set(gca, 'Xtick',[],'Ytick',[],'Color','w')
axis([-1 1 -1 1]);
axis equal
axis off
