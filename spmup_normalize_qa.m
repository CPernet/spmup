function spmup_normalize_qa(flags,Normalized,T1)

% implement different quality control 'tools' for
% fMRI spatially normalized data using SPM
%
% FORMAT normalize_qa ( = do it all)
%        normalize_qa(P,flags)
%
% INPUT flags is a structure listing the options available for qa
%             flags.average = 'on'; 
%                   --> create the average of normalized images
%                   --> if flags.T1 is on, also plot the outlines of the
%                   difference with the normalized T1
%             flags.T1 = 'on';
%                   --> compute the difference between the T1 (assumed
%                       normalized) and the normalized images - plot the outlines
%             flags.movie = 'on';
%             flags.AC = [x,y,z]; 
%                   --> generate three movies along the x, y, z plabes defined based on the AC point
%                   --> if flags.AC is empty, the user is asked to input this coordinate
%       Normalized is a matrix with the full name of files (see spm_select)
%       T1 is a vector with the full name of the normalized T1 file
%
% see also spm_select 
%
%  Cyril Pernet January 2013
% --------------------------------------------------------------------------
% Copyright (c) SPM U+ toolbox


spm('Defaults','fmri')
current = pwd;

% check flags
def_flags.T1 = 'on';
def_flags.average = 'on';
def_flags.movie = 'on';
if nargin == 0;
    flags = def_flags;
end

% update flags fields
if ~isfield(flags,'T1'); flags.T1 = 'off'; end
if ~isfield(flags,'average'); flags.average = 'off'; end
if ~isfield(flags,'movie'); flags.movie = 'off'; end

% check images as input
if nargin < 2
    if strcmp(flags.average,'on') 
        Normalized = spm_select(Inf,'image','select normalized images');
         V = spm_vol(Normalized); 
         if spm_check_orientations(V) ~= 1; error('Volumes are not all oriented the same way'); end
         if ~strcmp(V(1).descrip,'spm - 3D normalized'); warning('Volumes appear not to be SPM normalized ones'); end
    end
    if strcmp(flags.T1,'on')
        T1 = spm_select(1,'image','select normalized T1 image');
        T = spm_vol(T1); if ~strcmp(T.descrip,'spm - 3D normalized');
            warning('The T1 volume appears not to be SPM normalized'); 
        end
    end
elseif nargin < 3
    V = spm_vol(Normalized);
    if spm_check_orientations(V) ~= 1; error('Volumes are not all oriented the same way'); end
    if sum(findstr(V(1).descrip,'normalized'))==0 && sum(findstr(V(1).descrip,'normed'))==0
        warning('Volumes appear not to be SPM normalized ones'); end
    if strcmp(flags.T1,'on')
        T1 = spm_select(1,'image','select normalized T1 image');
        T = spm_vol(T1); if ~strncmp(T.descrip,'spm - 3D normalized',19);
            warning('The T1 volume appears not to be SPM normalized');
        end
    end
elseif nargin == 3
    V = spm_vol(Normalized);
    if spm_check_orientations(V) ~= 1; error('Volumes are not all oriented the same way'); end
    if sum(findstr(V(1).descrip,'normalized'))==0 && sum(findstr(V(1).descrip,'normed'))==0
        warning('Volumes appear not to be SPM normalized ones'); end
    T = spm_vol(T1);
    if sum(findstr(T.descrip,'normalized'))==0 && sum(findstr(T.descrip,'normed'))==0 
        warning('The T1 volume appears not to be SPM normalized'); end
end


[folder,name,ext] = fileparts(Normalized(1,:));
cd(folder); 

%% make average file
if strcmp(flags.average, 'on')
     disp('creating average normalized image')
     Images = spm_read_vols(V);
     Average = mean(Images,4);
     [folder,name,ext] = fileparts(V(1).fname);
     V(1).fname = [folder filesep 'normalized_mean_image' ext];
     spm_write_vol(V(1),Average);
end


%% outline between T1 and T2*
if strcmp(flags.T1, 'on')
    Ref = spm_read_vols(T);
    if strcmp(flags.average, 'off')
        Average = mean(spm_read_vols(V),4);
    end
    search = squeeze(Average(round(T.dim(1)/2),round(T.dim(2)/2),:));
    slices = find(search>10-3); slices = [slices(4):6:slices(end)];
    
    figure('Name','Average T1 vs T2*'); colormap('gray')
    for s=1:8
        A = flipud(squeeze(Ref(:,:,slices(s)))'); subplot(3,8,s); imagesc(A);
        if s ==4; title('Normalized T1 weighted image','Fontsize',14); end
        B = flipud(squeeze(Average(:,:,slices(s)))'); subplot(3,8,s+8); imagesc(B);
        if s ==4; title('Normalized averaged T2* weighted image','Fontsize',14); end
        C = abs((A>10-3) - (B>10-3)); subplot(3,8,s+16); imagesc(C);
        if s ==4; title('Difference in positive (mask) values','Fontsize',14); end
    end
    saveas(gcf,'outlines.fig','fig'); close(gcf)

end

%% if movie we need AC
if strcmp(flags.movie,'on')
    if ~isfield(flags,'AC')
        % display
        try
            Ref = spm_read_vols(T);
            spm_image('init',T)
            flags.AC = spm_input('enter a AC voxel coordinates',1);
        catch
            spm_image('init',V(1))
            flags.AC = spm_input('enter a AC voxel coordinates',1);
        end
        close all
    end
    
    for plane = 1:3
        figure; colormap('gray')
        X = squeeze(Images(flags.AC(1),:,:,:));
        Y = squeeze(Images(:,flags.AC(2),:,:));
        Z = squeeze(Images(:,:,flags.AC(3),:));

        % make video object
        if plane == 1
            vidObj = VideoWriter('sagital_normalized.avi');
        elseif plane == 2
            vidObj = VideoWriter('coronal_normalized.avi');
        elseif plane == 3
            vidObj = VideoWriter('axial_normalized.avi');
        end
        open(vidObj);
        set(gca,'NextPlot','replacechildren')
        
        % make pictures and get frames
        for i=1:size(Images,4)
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
disp('normalization QA done');
