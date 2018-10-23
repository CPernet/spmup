function new_file = spmup_realign_qa(P,flags)

% implement different quality control 'tools' for
% fMRI realigned data using SPM
%
% FORMAT realign_qa ( = do it all)
%        realign_qa(P,flags)
%
% INPUT flags is a structure listing the options available for qa
%             flags.motion_parameters = 'on'; 
%                   --> plot motion paramters and 1st derivatives
%                   --> generate new regressor.txt file with motion parameters and outliers detected using Carling's modified boxplot rule
%             flags.globals = 'on';
%                   --> compute and plot globals (see spm_globals)
%                   --> append the regressor.txt file with global outliers detected using Carling's modified boxplot rule
%             flags.volume_distance = 'on';
%                   --> compute the distance of each volume to the mean
%             flags.movie = 'on';
%             flags.AC = [x,y,z]; 
%                   --> generate three movies along the x, y, z planes defined based on the AC point
%                   --> if flags.AC is empty, the user is asked to input this coordinate
%       P is a matrix with the full name of files (see spm_select)
%
% see also spm_select spm_globals
%
% Cyril Pernet January 2013
% v2 added the new_files to work with the batch
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2015

%% validate inputs
spm('Defaults','fmri')
current = pwd;
 
% check flags
def_flags.motion_parameters = 'on';
def_flags.globals = 'on';
def_flags.volume_distance = 'on';
def_flags.movie = 'off';
if nargin == 0;
    flags = def_flags;
end

% update flags fields
if ~isfield(flags,'motion_parameters'); flags.motion_parameters = 'off'; end
if ~isfield(flags,'globals'); flags.globals = 'off'; end
if ~isfield(flags,'volume_distance'); flags.volume_distance = 'off'; end
if ~isfield(flags,'movie'); flags.movie = 'off'; end

% similar fields coming from the batch
if isfield(flags,'realign_QA_motion')
    if flags.realign_QA_motion == 1; flags.motion_parameters = 'on'; else flags.motion_parameters = 'off'; end
end
    
if isfield(flags,'realign_QA_globals')
    if flags.realign_QA_globals == 1; flags.globals = 'on'; else flags.globals = 'off'; end
end

if isfield(flags,'realign_QA_distance')
    if flags.realign_QA_distance ==1; flags.volume_distance = 'on'; else flags.volume_distance = 'off'; end
end

if isfield(flags,'realign_QA_movie')
    if flags.realign_QA_movie ==1; flags.movie = 'on'; end
end


% check P
if nargin <= 1
    % if globals, or distance, or movies we need P
    if strcmp(flags.globals,'on') || strcmp(flags.volume_distance,'on') || strcmp(flags.movie,'on')
        [P,check] = spm_select(Inf,'image','select realigned images');
        if check == 0
            return
        end
    end
end

if iscell(P)
    try
        P = cell2mat(P);
        V = spm_vol(P)
    catch
        for v=1:size(P,1)
            V(v) =spm_vol(P{v});
        end
    end
else
    V = spm_vol(P);
end
[folder,name,~] = fileparts(V(1).fname);
if strcmp(V(1).descrip, 'Warped') %in case we are dealing with a normalized image
    average = spm_select('FPList',folder,['^wmean' name(2:end)]);
else
    average = spm_select('FPList',folder,['^mean' name]);
end
cd(folder);

% if distance we need the mean image
if strcmp(flags.volume_distance,'on')
    if isempty(average)
        disp('mean image not found - computing one .. ')
        average = mean(spm_read_vols(V),4);
    else
        average =  spm_read_vols(spm_vol(average));
    end
end

% if movie we need AC
if strcmp(flags.movie,'on')
    if ~isfield(flags,'AC') && isempty(average)
        disp('mean image not found - computing one .. ')
        average = mean(spm_read_vols(V),4);
    end
    flags.AC = [size(average,1)/2, size(average,2)/2, size(average,3)/2];
end

 

%% look at motion parameters
moutlier_matrix = [];

if strcmp(flags.motion_parameters,'on')
    
    motion_file = dir('rp*.txt'); % SPM
    if isempty(motion_file)
        motion_file = dir('*_mcf.par'); % FSL
    end
    motion_param = load(motion_file.name);
    [n,p]=size(motion_param);
    
    disp('computing motion displacement for outliers ... ')
    % find outliers for high motion
    temp_motion = motion_param;
    temp_motion(:,4:6)= temp_motion(:,4:6)*180/pi; % use deg rather than radian
    derivatives = diff(temp_motion);
    
    delta = zeros(n,1);  % Mean square displacement in two scans
    for i = 2:n
        delta(i) = sqrt((temp_motion(i-1,1) - temp_motion(i,1))^2 + ...
            (temp_motion(i-1,2) - temp_motion(i,2))^2 +...
            (temp_motion(i-1,3) - temp_motion(i,3))^2 +...
            1.28*(temp_motion(i-1,4) - temp_motion(i,4))^2 +...
            1.28*(temp_motion(i-1,5) - temp_motion(i,5))^2 +...
            1.28*(temp_motion(i-1,6) - temp_motion(i,6))^2);
    end
    
    % robust outliers
    m_outliers = spm_comp_robust_outliers(delta);
    
    if sum(m_outliers) > 0
        moutlier_matrix = zeros(size(motion_param,1),sum(m_outliers));
        indices = find(m_outliers);
        for i=1:sum(m_outliers)
            moutlier_matrix(indices(i),i) = 1;
        end
    end
    
    % figure
    figure('Name','Motion outlier detection','Visible','On');
    set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized', 'outerposition',[0 0 1 1]); 
    subplot(3,2,1); plot(temp_motion(:,1:3),'LineWidth',3);
    axis tight; grid on; title('Translation','Fontsize',14)
    xlabel('scans'); ylabel('mm');
    subplot(3,2,2); plot(temp_motion(:,4:6),'LineWidth',3);
    axis tight; grid on; title('Rotation','Fontsize',14)
    xlabel('scans'); ylabel('degrees');
    subplot(3,2,3); plot(derivatives(:,1:3),'LineWidth',3);
    axis tight; grid on; title('scan to scan translation','Fontsize',14)
    xlabel('scans'); ylabel('mm');
    subplot(3,2,4); plot(derivatives(:,4:6),'LineWidth',3);
    axis tight; grid on; title('scan to scan rotation','Fontsize',14)
    xlabel('scans'); ylabel('degrees');
    subplot(3,2,5:6); plot(delta,'LineWidth',3); hold on;
    plot(m_outliers.*max(delta),'or','LineWidth',3);
    axis tight; grid on; title('Mean square displacement','FontSize',14);
    xlabel('scans'); ylabel('displacement');
    try
        saveas(gcf, 'motion outlier detection.fig','fig'); close(gcf)
    end

end


%% look at globals
goutlier_matrix = [];

if strcmp(flags.globals,'on')
    
    % find outliers in global mean intensity
    disp('computing globals for outliers ... ')
    glo = zeros(size(P,1),1);
    for s=1:size(V,2)
        glo(s) = spm_global(V(s));
    end
    glo = detrend(glo); % since in spm the data are detrended
    
    % robust outliers
    g_outliers = spm_comp_robust_outliers(glo);
    
    if sum(g_outliers) > 0
        goutlier_matrix = zeros(size(motion_param,1),sum(g_outliers));
        indices = find(g_outliers);
        for i=1:sum(g_outliers)
            goutlier_matrix(indices(i),i) = 1;
        end
    end
    
    % figure
    figure('Name','Globals outlier detectio','Visible','On');
    set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized', 'outerposition',[0 0 1 1]); 
    plot(glo,'LineWidth',3); hold on; plot(g_outliers.*(max(glo)-mean(glo))+mean(glo),'or','LineWidth',2);
    grid on; axis tight; title('Global intensity','FontSize',14); xlabel('scans'); ylabel('mean intensity')
    try
        saveas(gcf, 'globals outlier detection.fig','fig'); close(gcf)
    end

end


%% check distances to the mean and volume to volume
 dist_outliers_matrix = [];

if strcmp(flags.volume_distance, 'on')
    
    n = numel(V);
    distance_to_mean = zeros(n,1);
    distance_between = zeros(n-1,1);
    % mean_data = spm_read_vols(spm_vol(average));
    mask = average > 0;
    
    for i=1:n
        fprintf('computing mean square distances, volume %g/%g \n',i,n)
        vol = spm_read_vols(V(i));
        dist = (average(mask) - vol(mask)).^2;
        distance_to_mean(i) = sum(dist(~isnan(dist)))/numel(~isnan(dist));
        if i<n
            vol2 = spm_read_vols(V(i+1));
            dist = (vol2(mask) - vol(mask)).^2;
            distance_between(i) = sum(dist(~isnan(dist)))/numel(~isnan(dist));
        end
    end
    
    % robust outliers
    dist_outliers = spm_comp_robust_outliers(distance_between);
    
    figure('Name','Mean square distances between volumes','Visible','On');
    set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized', 'outerposition',[0 0 1 1]); 
    subplot(2,1,1); plot([1:n],distance_to_mean,'o','LineWidth',3);
    hold on; plot([1:n],distance_to_mean,'r','LineWidth',3);
    axis tight; grid on; xlabel('scans'); ylabel('square difference');
    title('Mean square distance to the average','FontSize',14);
    subplot(2,1,2); plot([1:n-1],[distance_between],'o','LineWidth',3);
    hold on;  plot([distance_between],'r','LineWidth',2);
    plot(dist_outliers.*(max(distance_between)-mean(distance_between))+mean(distance_between),'or','LineWidth',2);
    axis tight; grid on; xlabel('scans'); ylabel('square difference');
    title('Volume to volume mean square distance','FontSize',14);
    try
        saveas(gcf,'distances.fig','fig'); close(gcf)
    end
    
    if sum(dist_outliers) > 0
        dist_outliers_matrix = zeros(size(motion_param,1),sum(dist_outliers));
        dist_outliers = [0;dist_outliers]; indices = find(dist_outliers);
        for i=1:sum(dist_outliers)
            dist_outliers_matrix(indices(i),i) = 1;
        end
    end
end

%% update the motion parameter file

out = [moutlier_matrix goutlier_matrix dist_outliers_matrix]; 
if ~isempty(out)
    % check a volume is not indexed twice
    for j=1:size(out,2); w(j) = find(out(:,j)==1); end
    [~,index]=unique(w); out = out(:,index);
    multiple_regressors = [motion_param out];
else
    multiple_regressors = [motion_param];
end
save multiple_regressors.txt multiple_regressors -ascii

if nargout
    new_file = [pwd filesep 'multiple_regressors.txt'];
end




%% make movies
if strcmp (flags.movie, 'on')

    Images = spm_read_vols(V);
    if isnan(AC)
       AC(1) = round(size(Images,1) / 2);
       AC(2) = round(size(Images,2) / 2);
       AC(3) = round(size(Images,3) / 2);
    end
    X = squeeze(Images(flags.AC(1),:,:,:));
    Y = squeeze(Images(:,flags.AC(2),:,:));
    Z = squeeze(Images(:,:,flags.AC(3),:));
    clear Images
    
    for plane = 1:3
        figure; colormap('gray')
        
        % make video object
        if plane == 1
            vidObj = VideoWriter('sagital.avi');
        elseif plane == 2
            vidObj = VideoWriter('coronal.avi');
        elseif plane == 3
            vidObj = VideoWriter('axial.avi');
        end
        open(vidObj);
        set(gca,'NextPlot','replacechildren')
        
        % make pictures and get frames
        for i=1:n
            if plane == 1
                imagesc(X(:,:,i)');axis tight
            elseif plane == 2
                imagesc(Y(:,:,i)');axis tight
            elseif plane == 3
                imagesc(Z(:,:,i)');axis tight
            end
            title(['volume ' num2str(i)]);
            currFrame = getframe;
            writeVideo(vidObj,currFrame);
        end
        close(vidObj);close
    end

end

cd(current)
disp('realign QA done');

