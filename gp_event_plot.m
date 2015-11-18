function Y = gp_event_plot(varargin)

% plotting tool for group level event related response
% 
% FORMAT gp_event_plot
%               --> called via GUI, use the workspace to get varargin
%
%        Y = gp_event_plot(Coordinate,GpSPM,flag)
%               --> manual call allwing to obtain sets of responses
%
% INPUT Coordinate a 3*n set of voxel coordinates in MNI space
%       GpSPM the full (with path) name of the group analysis SPM.mat
%       flag indicates the coordinate space 'mm' or 'voxel', if not
%       specified user is prompted
%
% OUTPUT Y a structure with parameters information
%        Y.individual_responses = the fitted data per subject;
%        Y.individual_parameters = the coefficient (betas/cons) per subject;
%        Y.coordinate = the coordinate(s) used
%        Y.average.condition{n}.name = condition name(s);
%        Y.average.condition{n}.response = average response over subjects
%        Y.average.condition{n}.CI = 95% boostrap CI of the average reponse
%        
%        If no outpout Y is returned in ther workshaoce as gp_event and 
%        plots response(s) in a new window
%
% Cyril Pernet 13 Dec 2013.
% also now plot fitted data if boosted images are used and further
% plit the results per condition - Oct 2015
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
    cd(xSPM.swd);
    
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
disp('collecting data for all subjects ... '); clear tmp
for i=1:size(V,1) % for each subject
    
    % get the SPM.mat of this image
    [spath,name,ext]=fileparts(GpSPM.xY.P{i});
    try
        cd(spath);
        try
            load SPM
        catch SPM_not_located
            cd .. % if data in the dir above
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
        
        if strncmp(name,'boost',5)
            % compute the overall fitted response
            if strcmp(SPM.xBF.name,'hrf (with time derivative)') 
                tmp(1) = spm_get_data([SPM.swd filesep name(7:end) ext],Coordinate(:,c));
                add_one = num2str(eval(name(end-3:end)) + 1); newname = [name(7:end-length(add_one)) add_one];
                tmp(2) = spm_get_data([SPM.swd filesep newname ext],Coordinate(:,c));
                response{i,c} = SPM.xBF.bf(:,[1 2])*tmp';
            else % necesarilly time and dispersion
                tmp(1) = spm_get_data([SPM.swd filesep name(7:end) ext],Coordinate(:,c));
                add_one = num2str(eval(name(end-3:end)) + 1); newname = [name(7:end-length(add_one)) add_one];
                tmp(2) = spm_get_data([SPM.swd filesep newname ext],Coordinate(:,c));
                add_two = num2str(eval(name(end-3:end)) + 2); newname = [name(7:end-length(add_two)) add_two];
                tmp(3) = spm_get_data([SPM.swd filesep newname ext],Coordinate(:,c));
                response{i,c} = SPM.xBF.bf(:,[1 2 3])*tmp';
            end
        else
            % from here get the event related response
            response{i,c} = SPM.xBF.bf(:,1)*coef{i,c};
        end
        
        if c==1
            times = [0:SPM.xBF.dt:(length(SPM.xBF.bf)-1)*SPM.xBF.dt];
        else
            checktime = [0:SPM.xBF.dt:(length(SPM.xBF.bf)-1)*SPM.xBF.dt];
            if sum(checktime ~= times) ~=0
                times = [];
            end
        end
    end
end
clear Y tmp
cd(current);

%% outputs

disp('computing mean response and bootstrap 95% CI')
Y.individual_responses = response;
Y.individual_parameters = coef;
Y.coordinate = Coordinate;

% compute the mean response per condition
load SPM; index = 1;
for n=1:size(SPM.xX.name)
    if ~strncmp(SPM.xX.name{n},'subject',7)
        cname{index} = SPM.xX.name{n};
        index = index+1;
    end
end
Ncond = index-1;

clear tmp;
for c=1:size(Coordinate,2)
    % get all responses peer subjects
    for i=1:size(V,1)
        vv(i,:) = Y.individual_responses{i,c}; 
    end
    
    % get all response per condition and average
    for s =1:(size(V,1)/Ncond)
        for i=1:Ncond
            tmp{i}(:,:,s) = response{(s+i-1),c};
        end
    end
end

% now average across coordinates and subjects
for n=1:Ncond
    data = tmp{n}; data = squeeze(mean(data,2)); % average coordinates
    Y.average.condition{n}.name = cname{n};
    Y.average.condition{n}.response = mean(data,2); % average subjects
    
    % bootstrap avweraging the other way around to keep voxel variance
    data = tmp{n}; boot_data = NaN(size(data,1), size(data,2), 599);
    for b=1:599
        boot_data(:,:,b) = mean(data(:,:,randi(size(data,3),size(data,3),1)),3); % average over subjects
    end
    boot_data = sort(squeeze(mean(boot_data,2)),2); % average voxels and sort resamples
    Y.average.condition{n}.CI = [boot_data(:,15) boot_data(:,584)];
end

if nargout == 0
    assignin('base','gp_event',Y);
end

%% plot
if nargout == 0 
    figure('Name','Gp level evoked response','units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'Color','w','InvertHardCopy','off'); colormap('gray'); hold on;
    if isempty(times)
        times = [1:size(Y.average.condition{1}.response)];
    end
    mycolors = jet; mycolors = mycolors(1:64/Ncond:end,:);
    for n=1:Ncond
        plot(times,Y.average.condition{n}.response,'LineWidth',3,'Color',mycolors(n,:));
        plot(times,Y.average.condition{n}.CI(:,1),'LineWidth',2,'Color',mycolors(n,:));
        plot(times,Y.average.condition{n}.CI(:,2),'LineWidth',2,'Color',mycolors(n,:));
        fillhandle(n) = patch([times fliplr(times)], [Y.average.condition{n}.CI(:,1)' fliplr(Y.average.condition{n}.CI(:,2)')], mycolors(n,:));
        set(fillhandle(n),'EdgeColor',[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.8);
    end
    grid on; xlabel('PST (sec)','FontSize',14); 
    ylabel('Evoked Response (A.U.)','FontSize',14);
    if Ncond == 1
        title(sprintf('Average response\n at coord [%g %g %g]=%g [%g %g]',xyz(1),xyz(2),xyz(3)),'FontSize',16);
    else
        title(sprintf('Average responses per conditions \n at coord [%g %g %g]',xyz(1),xyz(2),xyz(3)),'FontSize',16);
    end
    legend(fillhandle,cname,'Location','SouthEast'); axis tight

    spm_results_ui('Clear',Fgraph);
    figure(Fgraph);
    subplot(2,2,4);
    for i=1:size(V,1); v(i) = coef{i,1}; end
    sv = sort(mean(v(randi(size(V,1),size(V,1),599))));
    L = mean(v)-sv(15); H = sv(584)-mean(v);
    bar(1,mean(v)); hold on; errorbar(1,mean(v),L,H,'LineWidth',3);
    axis([0.5 1.5 sv(15)-10/100*sv(15) sv(584)+10/100*sv(584)])
    title(sprintf('Average event related value and 95%% CI \n %g [%g %g]',mean(v), sv(15),sv(584)));

    subplot(2,2,3); hold on
    for n=1:Ncond
        plot(times,Y.average.condition{n}.response,'LineWidth',3,'Color',mycolors(n,:));
    end
    grid on; xlabel('PST'); ylabel('Event related resp.'); axis tight
    if Ncond == 1
        title(['Average response'],'FontSize',10);
    else
        title(['Average responses per conditions'],'FontSize',10);
    end
    
end
