function Y = gp_event_plot(varargin)

% plotting tool for group level event related response
% returns the adjusted response (ie the gp level adjusted beta+BF) as well
% as the average fitted data (ie the mean on indivudual modelled data)
%
% Note the if boosted data are used it recomputed the modelled data using
% the 3 basis functions while the adjusted data (from which you got stats)
% show the response wit the amplitude adjusted and the mean time 2 peak
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
%        Y.individual_parameters = the coeffient (betas) per subject;
%        Y.individual_adjusted_parameters = the coeffient adjusted by the
%                                    model (this is where stats come from)
%        Y.individual_responses = the fitted data per subejct;
%        Y.individual_adjusted_responses = the modelled data using adjusted coef;
%        Y.coordinate = the coordinate(s) used
%        Y.average.condition{n}.name = condition name(s);
%        Y.average.condition{n}.response = average response over subjects
%        Y.average.condition{n}.CI = 95% boostrap CI of the average reponse
%        if boosted data are used Y.individual_estimated_time_to_peak (and
%        the adjusted response if based in this)
%        if con imagaes are used Y.individual_beta_coef = beta coef combined in the contrast;
%
%        If no outpout if specificied (ie GUI), it also plots the response
%        in a new window.
%
% Cyril Pernet January 2016

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
        % vox=xyz'*tmp(1:3,1:3) + tmp(1:3,4)';
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
GpSPM = SPM;
V = spm_vol(GpSPM.xY.P);


%-Get coef with whittening (ie where we got the stats from)
%-----------------------------------------------------------
disp('collecting data for all subjects ... '); clear tmp

% the coef of the input images are:
y = spm_get_data(V,Coordinate);
if sum(y) == 0
    error('all estimates are 0 ???, typically an error of input mm vs voxel')
end

if any(isnan(y))
    error('at least some estimates are NaN ?? check input coordinate')
end

% the adjusted coef of the input images at the gp level are:
yy = GpSPM.xX.W*y;


% reconstruct the response
% -------------------------

% hrf parameters
p(1) = 6;
p(2) = 16;
p(3) = 1;
p(4) = 1;
p(5) = 6;
p(6) = 0;
p(7) = 32;

for i=1:size(V,1) % for each image/subject
    % get file info
    [spath,name,ext]=fileparts(GpSPM.xY.P{i});
    
    % get the SPM.mat of this image
    cd(fileparts(GpSPM.xY.P{i}));
    try
        try
            load SPM; 
        catch
            cd ..
            load SPM
        end
    catch
        error('Can''t locate SPM.mat');
    end
    
    % check SPM.mat path
    if exist([SPM.swd filesep 'SPM.mat'],'file')
        SPMPath = SPM.swd;
    else
        SPMPath = pwd;
        if i==1
            warndlg('couldn''t locate SPM.mat based on saved info, using local directory structure','SPM.swd error')
        end
        fprintf('the SPM.mat of image %g contains invalid path information, this often happens when moving data\n', i)
    end

            
    % compute for each coordinate
    for c=1:size(Coordinate,2)
        
        % the coef of the input image and event model + time adjustment
        coef{i,c} = y(i,:); adjusted_coef{i,c} = yy(i,:);
        
        if strncmp(name,'boost',5) || strncmp(name,'sboost',6) % boosted parameter estimates
            if sum(findstr(name,'beta')) ~=0
                p(1) = spm_get_data(spm_vol([spath filesep 'T2P' num2str(eval(name(end-3:end))) ext]),Coordinate(:,c));
            elseif sum(findstr(name,'con')) ~=0
                beta_indices = find(SPM.xCon(eval(name(end-3:end))).c); % which betas were combined in the constrast
                for comb = 1:length(beta_indices)
                    t2p(comb) = spm_get_data(spm_vol([spath filesep 'T2P' num2str(beta_indices(comb)) ext]),Coordinate(:,c));
                end
                p(1) = mean(t2p);
            end
            hrf = spm_hrf(SPM.xBF.dt,p,SPM.xBF.T);
            adjusted_response{i,c} = hrf*adjusted_coef{i,c};
            estimated_time_to_peak{i,c} = p(1);
        else % standard parameter estimates
            adjusted_response{i,c} = SPM.xBF.bf(:,1)*adjusted_coef{i,c};
        end
        
        
        % the timing is
        if c==1
            times = [0:SPM.xBF.dt:(length(SPM.xBF.bf)-1)*SPM.xBF.dt];
        else
            checktime = [0:SPM.xBF.dt:(length(SPM.xBF.bf)-1)*SPM.xBF.dt];
            if sum(checktime ~= times) ~=0
                times = [];
            end
        end
        
        % the fitted responses per subject (non adjusted)
        if strncmp(name,'boost',5) || strncmp(name,'sboost',6)
            
            if strcmp(SPM.xBF.name,'hrf (with time derivative)')
                if sum(findstr(name,'beta')) ~=0
                    try
                        tmp(1) = spm_get_data([SPMPath filesep name(7:end) ext],Coordinate(:,c)); % get the coresponding hrf value
                        add_one = num2str(eval(name(end-3:end)) + 1); newname = [name(7:end-length(add_one)) add_one]; % get the time derivative beta
                        tmp(2) = spm_get_data([SPMPath filesep newname ext],Coordinate(:,c));
                    catch
                        tmp(1) = spm_get_data([SPMPath filesep name(8:end) ext],Coordinate(:,c)); 
                        add_one = num2str(eval(name(end-3:end)) + 1); newname = [name(8:end-length(add_one)) add_one]; 
                        tmp(2) = spm_get_data([SPMPath filesep newname ext],Coordinate(:,c));
                    end
                    response{i,c} = SPM.xBF.bf(:,[1 2])*tmp';
                elseif sum(findstr(name,'con')) ~=0
                    beta_indices = find(SPM.xCon(eval(name(end-3:end))).c); % which betas were combined in the constrast
                    for comb = 1:length(beta_indices)
                        if beta_indices(comb)<10
                            betaname = ['beta_000' num2str(beta_indices(comb))];
                        elseif beta_indices(comb)<100
                            betaname = ['beta_00' num2str(beta_indices(comb))];
                        elseif beta_indices(comb)<1000
                            betaname = ['beta_0' num2str(beta_indices(comb))];
                        else
                            betaname = ['beta_' num2str(beta_indices(comb))];
                        end
                        tmp(comb,1) = spm_get_data([SPMPath filesep betaname ext],Coordinate(:,c));
                        add_one = num2str(eval(betaname(end-3:end)) + 1); newname = [betaname(1:end-length(add_one)) add_one];
                        tmp(comb,2) = spm_get_data([SPMPath filesep newname ext],Coordinate(:,c));
                    end
                    beta_coef{i,c} = tmp; % keep orignal betas
                    tmp = mean(tmp,1); % average over beta value
                    response{i,c} = SPM.xBF.bf(:,[1 2])*tmp';
                end
                
            else % necesarilly time and dispersion
                if sum(findstr(name,'beta')) ~=0
                    try
                        tmp(1) = spm_get_data([SPMPath filesep name(7:end) ext],Coordinate(:,c));
                        add_one = num2str(eval(name(end-3:end)) + 1); newname = [name(7:end-length(add_one)) add_one];
                        tmp(2) = spm_get_data([SPMPath filesep newname ext],Coordinate(:,c));
                        add_two = num2str(eval(name(end-3:end)) + 2); newname = [name(7:end-length(add_two)) add_two]; % get the dispersion derivative beta
                        tmp(3) = spm_get_data([SPMPath filesep newname ext],Coordinate(:,c));
                    catch
                        tmp(1) = spm_get_data([SPMPath filesep name(8:end) ext],Coordinate(:,c));
                        add_one = num2str(eval(name(end-3:end)) + 1); newname = [name(8:end-length(add_one)) add_one];
                        tmp(2) = spm_get_data([SPMPath filesep newname ext],Coordinate(:,c));
                        add_two = num2str(eval(name(end-3:end)) + 2); newname = [name(8:end-length(add_two)) add_two];
                        tmp(3) = spm_get_data([SPMPath filesep newname ext],Coordinate(:,c));
                    end
                    response{i,c} = SPM.xBF.bf(:,[1 2 3])*tmp';
                elseif sum(findstr(name,'con')) ~=0
                    beta_indices = find(SPM.xCon(eval(name(end-3:end))).c);
                    for comb = 1:length(beta_indices)
                        if beta_indices(comb)<10
                            betaname = ['beta_000' num2str(beta_indices(comb))];
                        elseif beta_indices(comb)<100
                            betaname = ['beta_00' num2str(beta_indices(comb))];
                        elseif beta_indices(comb)<1000
                            betaname = ['beta_0' num2str(beta_indices(comb))];
                        else
                            betaname = ['beta_' num2str(beta_indices(comb))];
                        end
                        tmp(comb,1) = spm_get_data([SPMPath filesep betaname ext],Coordinate(:,c));
                        add_one = num2str(eval(betaname(end-3:end)) + 1); newname = [betaname(1:end-length(add_one)) add_one];
                        tmp(comb,2) = spm_get_data([SPMPath filesep newname ext],Coordinate(:,c));
                        add_two = num2str(eval(betaname(end-3:end)) + 2); newname = [betaname(1:end-length(add_two)) add_two];
                        tmp(comb,3) = spm_get_data([SPMPath filesep newname ext],Coordinate(:,c));
                    end
                    beta_coef{i,c} = tmp; % keep orignal betas
                    tmp = mean(tmp,1); % average over beta value
                    response{i,c} = SPM.xBF.bf(:,[1 2 3])*tmp';
                end
            end
            
        else
            % from here get the event related response
            response{i,c} = SPM.xBF.bf(:,1)*coef{i,c};
        end
        
    end
end
clear Y tmp
cd(current);

%% outputs

disp('computing mean response and bootstrap 95% CI')
Y.coordinate = Coordinate;
Y.individual_responses = response;
Y.individual_adjusted_responses = adjusted_response;
Y.individual_parameters = coef;
Y.individual_adjusted_parameters = adjusted_coef;
if exist('estimated_time_to_peak','var')
    Y.individual_estimated_time_to_peak = estimated_time_to_peak;
end
if exist('beta_coef','var')
    Y.individual_beta_coef = beta_coef;
end


% compute the mean response per condition
index = 1; repeated_measure = 'no';
for n=1:size(GpSPM.xX.name,1)
    if ~strncmpi(GpSPM.xX.name{n},'subject',7)
        cname{index} = GpSPM.xX.name{n};
        index = index+1;
    else
        repeated_measure = 'yes';
    end
end
Ncond = index-1;
N = size(Y.individual_responses,1)/Ncond;


% now average across coordinates 
clear tmp t2p
warning off
for n=1:Ncond
    if exist('estimated_time_to_peak','var')
        t2p = cell2mat(estimated_time_to_peak(find(GpSPM.xX.X(:,n)))');
        if numel(size(t2p)) == 3
            t2p  = squeeze(t2p(tmp,3));
        end
        Y.average.condition{n}.time_to_peak = mean(t2p,2);
    end
    
    
     % ---- 
    
    tmp = cell2mat(Y.individual_adjusted_parameters(find(GpSPM.xX.X(:,n)))');
    data = cell2mat(Y.individual_adjusted_responses(find(GpSPM.xX.X(:,n)))');
    if numel(size(data)) == 3
        tmp = squeeze(mean(tmp,3)); % average coordinates
        t2p  = squeeze(t2p(tmp,3));
        data = squeeze(mean(data,3)); 
    end
    Y.adjusted_average.condition{n}.name = cname{n};
    Y.adjusted_average.condition{n}.parameters = mean(tmp,2);% average subjects
    Y.adjusted_average.condition{n}.response = mean(data,2); 
    
    fprintf('computing average adjusted response boostrap CI condition %g \n',n)
    boot_data = NaN(size(data,1), 599); % bootstrap
    for b=1:599
        go = 0; 
        while go == 0
            tmp = randi(size(data,2),size(data,2),1);
            if length(unique(tmp)) > 3
                rand_table(:,b) = tmp; go = 1;
            end
        end
        boot_data(:,b) = squeeze(mean(data(:,rand_table(:,b)),2)); % resample subjects and average
    end
   
    boot_data=sort(boot_data,2);
    for t=1:size(boot_data,1)
        if Y.adjusted_average.condition{n}.response(t)<=boot_data(t,15)
            Y.adjusted_average.condition{n}.CI = [boot_data(:,15) boot_data(:,584)];
        else
            Y.adjusted_average.condition{n}.CI = [boot_data(:,584) boot_data(:,15)];
        end
    end
        
    % ---- 
    tmp = cell2mat(Y.individual_parameters(find(GpSPM.xX.X(:,n)))');
    data = cell2mat(Y.individual_responses(find(GpSPM.xX.X(:,n)))');
    if numel(size(data)) == 3
        tmp = squeeze(mean(tmp,3)); % average coordinates
        t2p  = squeeze(t2p(tmp,3));
        data = squeeze(mean(data,3)); 
    end
    Y.average.condition{n}.name = cname{n};
    Y.average.condition{n}.parameters = mean(tmp,2);% average subjects
    Y.average.condition{n}.response = mean(data,2); 
   
    fprintf('computing average response boostrap CI condition %g \n',n)
    boot_data = NaN(size(data,1), 599); % bootstrap
    for b=1:599
        boot_data(:,b) = squeeze(mean(data(:,rand_table(:,b)),2)); % resample subjects and average
    end
   
    boot_data=sort(boot_data,2);
    for t=1:size(boot_data,1)
        if Y.average.condition{n}.response(t)<=boot_data(t,15)
            Y.average.condition{n}.CI = [boot_data(:,15) boot_data(:,584)];
        else
            Y.average.condition{n}.CI = [boot_data(:,584) boot_data(:,15)];
        end
    end
end
warning on

Y.individual_parameters = reshape(cell2mat(coef),Ncond,N)';
Y.individual_adjusted_parameters = reshape(cell2mat(adjusted_coef),Ncond,N');


%% plot
      
    % average adjusted responses and plot 95% CI
    figure('Name','Gp level evoked response','units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'Color','w','InvertHardCopy','off'); colormap('gray'); hold on;
    if isempty(times)
        times = [1:size(Y.average.condition{1}.response)];
    end
    mycolors = rst_colour_maps(Ncond) ; % cubehelixmap('semi_continuous',Ncond+2);
    
    for n=1:Ncond
        plot(times,Y.average.condition{n}.response,'LineWidth',3,'Color',mycolors(n,:)); hold on
        plot(times,Y.average.condition{n}.CI(:,1)','LineWidth',0.5,'Color',mycolors(n,:));
        plot(times,Y.average.condition{n}.CI(:,2)','LineWidth',0.5,'Color',mycolors(n,:));
        fillhandle(n) = patch([times fliplr(times)], [Y.average.condition{n}.CI(:,1)' fliplr(Y.average.condition{n}.CI(:,2)')], mycolors(n,:));
        set(fillhandle(n),'EdgeColor',mycolors(n,:),'FaceAlpha',0.2,'EdgeAlpha',0.2);
    end
    grid on; xlabel('PST (sec)','FontSize',14);
    ylabel('Evoked Response (A.U.)','FontSize',14);
    if Ncond == 1
        title(sprintf('Average response\n at coord [%g %g %g]=%g [%g %g]',xyz(1),xyz(2),xyz(3)),'FontSize',16);
    else
        title(sprintf('Average responses per conditions \n at coord [%g %g %g]',xyz(1),xyz(2),xyz(3)),'FontSize',16);
    end
    legend(fillhandle,cname,'Location','SouthEast'); axis tight
    
    % get the individual parameters
    [mean_param,CI_param]=rst_data_plot(Y.individual_parameters,'estimator','mean','newfig','yes');
    title('Parameter estimates for each condition')
    for n=1:Ncond
        Y.average.condition{n}.param_CI = CI_param(:,n);
    end
    
if nargout == 0 % called via GUI

    % update the SPM figure too
    spm_results_ui('Clear',Fgraph);
    figure(Fgraph);
    subplot(2,2,4);
    L = mean(mean_param-mean(CI_param(1,:))); H = mean(mean_param-mean(CI_param(2,:)));
    bar(1,mean(mean_param)); hold on; errorbar(1,mean(mean_param),L,H,'LineWidth',3);
    title(sprintf('Average adjusted coef and 95%% CI \n %g [%g %g]',mean(mean_param), L, H));
    grid on; box on
    
    subplot(2,2,3); hold on
    for n=1:Ncond
        plot(times,Y.average.condition{n}.response,'LineWidth',3,'Color',mycolors(n,:));
    end
    grid on; xlabel('PST'); ylabel('Event related resp.'); axis tight
    if Ncond == 1
        title(['Average response'],'FontSize',10);
    else
        title(['Average responses'],'FontSize',10);
    end
    
    assignin('base','gp_event',Y);
end
