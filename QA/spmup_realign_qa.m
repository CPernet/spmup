function [new_files,FD,glo] = spmup_realign_qa(P,varargin)

% implement different quality control 'tools' for
% fMRI realigned data using SPM
%
% FORMAT spmup_realign_qa ( = do it all)
%        spmup_realign_qa(P,flags)
%
% INPUT P indicate the timeseries to load. Can be 
%         - path to a 4D nifti
%         - cellstring of path to 3D niftis
%         - a 4D array x, y, z, t of the timeseries data
%
%       Options are:
%
%       'Framewise displacement': 'on' (default) or 'off'
%              --> plots displacement (FD and RMS)                   
%              --> generates/appends new design.txt file with displacement outliers
%       'Voltera': 'off' (default) or 'on'
%              --> clompute the motion parameter expansion
%              --> generates/appends new design.txt file with displacement outliers
%       'Globals': 'on' (default) or 'off'
%              --> computes and plot globals (see spm_globals)
%              --> generates/appends new design.txt file with global outliers
%       'Movie': 'on' (default) or 'off'
%       'Coordinate': a coordinate, like [46 64 37];
%              --> generates three movies along the x, y, z planes defined based on the coordinate
%              --> if empty, takes the middle of the volume
%       'Average': where to get the mean image
%       'Figure': : 'on', 'save' (default - save on drive but not visible) or 'off'
%
% Cyril Pernet 
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2019

%% validate inputs
spm('Defaults','fmri')


% check options
FramewiseDisplacement = 'on';
Radius                = 50;
Voltera               = 'off';
Globals               = 'on'; 
Movie                 = 'on';
Coordinates           = [];
fig                   = 'save';

if nargin >1
   for v=1:(nargin-1)
      if strcmpi(varargin{v},'Radius')
           Radius                = varargin{v+1};
      elseif strcmpi(varargin{v},'Framewise Displacement')
           FramewiseDisplacement = varargin{v+1};
      elseif strcmpi(varargin{v},'Voltera')
           Voltera               = varargin{v+1};
      elseif strcmpi(varargin{v},'Globals')
           Globals               = varargin{v+1};
      elseif strcmpi(varargin{v},'Movie')
           Movie                 = varargin{v+1};
      elseif strcmpi(varargin{v},'Coordinate')
           Coordinates           = varargin{v+1};
      elseif strcmpi(varargin{v},'Figure')
           fig                   = varargin{v+1};
      end
   end
end

% check P
if nargin == 0
    [P,sts] = spm_select(Inf,'image' ,'Select your fMRI time series',{},pwd,'.*',Inf);
    if sts == 0
        return
    end
    V = spm_vol(P);
    if strcmp(Movie,'on')
        % bypass orientation check allowing to look at raw data
        N = numel(V);
        Y = zeros([V(1).dim(1:3),N]);
        for i=1:N
            for p=1:V(1).dim(3)
                Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
            end
        end
    end
else
    if ischar(P)
        V = spm_vol(P);
        if strcmp(Movie,'on')
            N = numel(V);
            Y = zeros([V(1).dim(1:3),N]);
            for i=1:N
                for p=1:V(1).dim(3)
                    Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
                end
            end
        end
    elseif iscell(P)
        for v=size(P,1):-1:1
            V(v) =spm_vol(P{v});
        end
        if strcmp(Movie,'on')
            for i=numel(V):-1:1
                for p=V(1).dim(3):-1:1
                    Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
                end
            end
        end
    else
        if numel(size(P)) == 4 % this is already data in
            Y = P;
            MotionParameters = 'off'; 
            FramewiseDisplacement = 'off';
            Voltera = 'off'; 
            Globals = 'on'; 
            Movie = 'on';
            fig = 'save';
        else
            error('input data are not char nor 4D data matrix, please check inputs')
        end
    end
end

% if movie we need AC
if strcmp(Movie,'on')
    if isempty(Coordinates)
        if exist('V','var')
            Coordinates = round(V(1).dim/2); % default
        else
            Coordinates = round(size(Y)/2);
            Coordinates = Coordinates(1:3);
        end
    end
end

findex              = 1; % index for the new_files variables
[filepath,filename] = fileparts(V(1).fname);
motion_file         = dir(fullfile(filepath,'rp*.txt'));
if size(motion_file,1)>1
    % remove eventual prefix in case we need to match 
    % the realignement parameters with a realigned file
    % that was prefixed
    unprefixed_filename = filename(3:end);
    motion_file = motion_file(arrayfun(@(x) contains(x.name, unprefixed_filename), motion_file)).name;
    motion_file = fullfile(filepath, motion_file);
else
    motion_file = fullfile(filepath, motion_file.name);
end

motion = load(motion_file);
if strcmpi(Voltera,'on')
    D      = diff(motion,1,1); % 1st order derivative
    D      = [zeros(1,6); D];  % set first row to 0
    motion = [motion D motion.^2 D.^2];
end

FD = [];
if strcmpi(FramewiseDisplacement,'on')
    FD = spmup_FD(motion_file, 'Radius', Radius, 'Figure', fig);
end

glo = [];
if strcmpi(Globals,'on')   
    glo = compute_globals(V);
    plot_globals(glo, fig);
end

%% make the design.txt file

design = motion; %#ok<*NASGU>
if strcmpi(FramewiseDisplacement,'on') || strcmpi(Globals,'on')
    % we do not censor the motion parameters 
    % but we keep them for the design
    data = [FD glo];
    design = spmup_censoring(data);
    design = [motion design];
end

save(fullfile(filepath,[filename '_design.txt']), 'design','-ascii')
new_files{findex} = fullfile(filepath,[filename '_design.txt']);

add_side_car(new_files{findex}, Voltera, FramewiseDisplacement, Globals);

findex = findex +1;

%% make movies
if strcmpi(Movie, 'on')
    new_files{findex} = spmup_movie(Y, 'coordinates', Coordinates,...
        'filename', fullfile(filepath,filename), 'showfig', 'off');
end

end

function add_side_car(new_file, Voltera, FramewiseDisplacement, Globals)

    opts.indent = ' '; % for spm_jsonwrite

    % save info about column headers
    options = struct('Voltera', Voltera, ...
                     'FramewiseDisplacement', FramewiseDisplacement, ...
                     'Globals', Globals);
    metadata.Columns = column_headers(options);
    
    all_regressors = spm_load(new_file);
    nb_censoring_regressors = size(all_regressors, 2) - numel(metadata.Columns);
    
    for i = 1:numel(1:nb_censoring_regressors)
        metadata.Columns{end+1} = sprintf('outlier_%04.0f', i);
    end
    spm_jsonwrite(spm_file(new_file, 'ext', '.json'), metadata, opts)
end

function headers = column_headers(options)
  
  headers = {'trans_x'; ...
            'trans_y'; ...
            'trans_z'; ...
            'rot_x'; ...
            'rot_y'; ...
            'rot_z'};
  
  if strcmpi(options.Voltera, 'on')
    
    headers = {'trans_x'; ...
              'trans_y'; ...
              'trans_z'; ...
              'rot_x'; ...
              'rot_y'; ...
              'rot_z'; ...
              'trans_x_derivative1'; ...
              'trans_y_derivative1'; ...
              'trans_z_derivative1'; ...
              'rot_x_derivative1'; ...
              'rot_y_derivative1'; ...
              'rot_z_derivative1'; ...
              'trans_x_power2'; ...
              'trans_y_power2'; ...
              'trans_z_power2'; ...
              'rot_x_power2'; ...
              'rot_y_power2'; ...
              'rot_z_power2'; ...
              'trans_x_derivative1_power2'; ...
              'trans_y_derivative1_power2'; ...
              'trans_z_derivative1_power2'; ...
              'rot_x_derivative1_power2'; ...
              'rot_y_derivative1_power2'; ...
              'rot_z_derivative1_power2'};
  end
  
  if strcmpi(options.FramewiseDisplacement, 'on')
       headers{end+1} = 'FramewiseDisplacement';
  end
  if strcmpi(options.Globals, 'on')
      headers{end+1} = 'Globals';
  end
    
end

function glo = compute_globals(V)
    glo = zeros(length(V),1);
    for s=1:length(V)
        glo(s) = spm_global(V(s));
    end
    glo        = spm_detrend(glo,1); % since in spm the data are detrended
end

function plot_globals(glo, fig)

        if ~strcmpi(fig,'off')
            return
        end
        
        figure('Name','Globals outlier detection')
        
        if strcmpi(fig,'on')
            set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
        else
            set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1], ...
                    'visible','off')
        end
        
        plot(glo,'LineWidth',3); 
        hold on;
        
        g_outliers = spmup_comp_robust_outliers(glo, 'Carling');
        tmp = g_outliers.*glo; tmp(tmp==0)=NaN;
        plot(tmp,'or','LineWidth',3); 
        
        grid on;
        axis tight; xlabel('scans');
        ylabel('mean intensity');
        title('Global intensity');
        
        if strcmpi(fig,'save')
            if exist(fullfile(filepath,'spmup_QC.ps'),'file')
                print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(filepath,'spmup_QC.ps'));
            else
                print (gcf,'-dpsc2', '-bestfit', fullfile(filepath,'spmup_QC.ps'));
            end
            close(gcf)
        end
end
