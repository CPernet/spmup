function Y = gp_event_plot(varargin)

% plotting tool for group level event related response
% 
% FORMAT gp_event_plot
%               --> called via GUI, use the workspace to get varargin
%
%        Y = gp_event_plot(Coordinate,GpSPM,flag)
%               --> manual call allwong to obtain sets of responses
%
% INPUT Coordinate a 3*n set of voxel coordinates in MNI space
%       GpSPM the full (with path) name of the group analysis SPM.mat
%       flag indicates the coordinate space 'mm' or 'voxel', if not
%       specified user is prompted
%
% OUTPUT Y a cell array of size n containing the group level event related 
%        response
%        
%        If no outpout if specificied (ie GUI), it also plots the response
%        in a new window.
%
% Cyril Pernet 13 Dec 2013.

%% check inputs

flag = [];
if nargin == 0
    
    Finter = spm_figure('GetWin','Interactive');
    Fgraph = spm_figure('GetWin','Graphics');
    hReg = evalin('base','hReg');
    xyz = [];
    xSPM = evalin('base','xSPM');
    GpPath = xSPM.swd;
    xyz = spm_XYZreg('NearestXYZ',...
        spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
    spm_XYZreg('SetCoords',xyz,hReg); % update GUI location
    
    if isempty(xyz)
        error('can''t find the coordinates')
    else
        tmp = inv(xSPM.Vspm.mat);
        Coordinate = tmp(1:3,:)*[xyz' 1]'; % change to voxel space
    end
    
elseif nargin ==2 || nargin ==3
    xyz = varargin{1};
    [GpPath name ext] = fileparts(varargin{2});
    
    if nargin ==3
        flag = varargin{3};
    else
        flag = spm_input('are coordinate in voxel?','!+1','b',{'mm','voxel'});
    end
    
    if strcmp(flag,'mm')
        cd(GpPath); load SPM
        tmp = inv(SPM.xVol.M);
        C = ones(4,size(xyz,2));
        C(1:3,:) = xyz;
        Coordinate = tmp(1:3,:)*C;
    elseif strcmp(flag,'voxel')
        Coordinate = xyz;
    else
        error('can''t find the coordinate space mm or voxel ????')
    end
else
    error('wrong input arguments')
end

%% get information needed to extract data and reconstruct the parametric response

% which images were used
current = pwd; cd(GpPath); load SPM
GpSPM = SPM; V = spm_vol(GpSPM.xY.P);

%-Get coef and reconstruct the response
%-----------------------------------------------------------------------
for i=1:size(V,1) % for each subject
    
    % get the SPM.mat of this image
    [spath,name,ext]=fileparts(GpSPM.xY.P{i});
    try
        cd(spath);
        try
            load SPM
        catch SPM_not_located
            cd .. % if data in a subdir
            load SPM
        end
    catch SPM_not_located
        [t,sts] = spm_select(1,'mat','SPM.mat not found',[],spath,[],[]);
        if sts == 1
            load(t)
        else
            return
        end
    end
    
    for c=1:size(Coordinate,2) % for each coordinate
        
        % the coef of the input image
        coef{i,c} = spm_get_data(V{i},Coordinate(:,c));
        
        % from here get the event related response
        response{i,c} = SPM.xBF.bf(:,1)*coef{i,c};
    end
end
clear Y tmp
cd(current);

%% outputs

Y.individual_responses = response;
Y.individual_parameters = coef;
Y.coordinate = Coordinate;

for c=1:size(Coordinate,2)
    for i=1:size(V,1)
        tmp(:,:,i) = response{i,c};
    end
    Y.average{c} = mean(tmp,3);
end

if nargout == 0
    assignin('base','gp_event',Y);
end

%% plot
if nargout == 0 
    spm_results_ui('Clear',Fgraph);
    figure(Fgraph);
    subplot(2,2,4);
    for i=1:size(V,1); v(i) = coef{i,1}; end
    sv = sort(mean(v(randi(size(V,1),size(V,1),599))));
    L = mean(v)-sv(15); H = sv(584)-mean(v);
    bar(1,mean(v)); hold on; errorbar(1,mean(v),L,H,'LineWidth',3);
    axis([0.5 1.5 0 sv(584)+10/100*sv(584)])
    title(sprintf('Average event related value and 95%% CI \n %g [%g %g]',mean(v), sv(15),sv(584)));

    subplot(2,2,3); plot(Y.average{1},'LineWidth',3); 
    for i=1:size(V,1); vv(i,:) = Y.individual_responses{i,1}; end
    clear tmp; for b=1:599
        tmp(:,:,b) = vv(randi(size(V,1),size(V,1),1),:);
    end
    svv = sort(squeeze(mean(tmp,1)),2); clear tmp; tmp = [1:size(svv,1)];
    hold on; plot(svv(:,15),'--','LineWidth',2); plot(svv(:,584),'--','LineWidth',2);
    grid on; xlabel('PST'); ylabel('Event related resp.');
    title(['Average response across subjects'],'FontSize',10);
    
    figure('Name','Gp level event related response')
    plot(Y.average{1},'LineWidth',3); hold on; plot(svv(:,15),'r','LineWidth',2); plot(svv(:,584),'r','LineWidth',2);
    fillhandle = patch([tmp fliplr(tmp)], [svv(:,15)' fliplr(svv(:,584)')], [1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.8);
    grid on; xlabel('PST'); ylabel('Event related resp.');
    title(sprintf('Average response across subjects \n at coord [%g %g %g]=%g [%g %g]',xyz(1),xyz(2),xyz(3),mean(v), sv(15),sv(584)),'FontSize',10);
end
