function jobs = spmup_first_level_qa(varargin)

% routine calling spmup_basics and spmup_realign_qa 
%
% INPUT jobs = spmup_first_level_qa
%       jobs = spmup_first_level_qa(time_series,options)
%
%       time_series is the fMRI set to analyze
%       options define various choices to be used:
%       'Motion Parameters': 'on' (default) or 'off'
%              --> plots motion parameters and 1st derivatives
%       'Radius': 57 (default)
%              --> the average head size to compute rotation displacement
%       'Framewise displacement': 'on' (default) or 'off'
%              --> plots displacement (FD and RMS)                   
%              --> generates/appends new design.txt file with displacement outliers
%       'Voltera': 'off' (default) or 'on'
%              --> compute the motion parameter expansion
%              --> generates/appends new design.txt file with motion derivatives and squares
%       'Globals': 'on' (default) or 'off'
%              --> computes and plot globals (see spm_globals)
%              --> generates/appends new design.txt file with global outliers
%       'Movie': 'on' (default) or 'off'
%       'Coordinate': a coordinate, like [46 64 37];
%              --> generates three movies along the x, y, z planes defined based on the coordinate
%       'Basics': 'on' (default) or 'off'
%              --> creates the average and std of images
%       'Figure': 'on', 'save' (default - save on drive but not visible) or 'off' 
% 
% OUTPUT jobs is a structure with the different jobs performed
%
% Cyril Pernet 
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2019

% defaults
spm('defaults','FMRI');

%% check inputs
if nargin == 0    
    [time_series,sts] = spm_select(Inf,'image' ,'Select your fMRI time series',{},pwd,'.*',Inf);
    if sts == 0
        return
    end
    
elseif nargin == 1
    time_series = varargin{1};
end


%% options
MotionParameters      = 'on';
Radius                = 57;
FramewiseDisplacement = 'on';
Voltera               = 'off';
Globals               = 'on'; 
Movie                 = 'on';
Basics                = 'on';
Figure                = 'save';

if nargin >1
   time_series = varargin{1};
   for v=1:nargin
      if strcmpi(varargin{v},'Motion Parameters')
           MotionParameters = varargin{v+1};
      elseif strcmpi(varargin{v},'Radius')
           Radius = varargin{v+1};
      elseif strcmpi(varargin{v},'Framewise Displacement')
           FramewiseDisplacement = varargin{v+1};
      elseif strcmpi(varargin{v},'Voltera')
           Voltera = varargin{v+1};
      elseif strcmpi(varargin{v},'Globals')
           Globals = varargin{v+1};
      elseif strcmpi(varargin{v},'Movie')
           Movie = varargin{v+1};
      elseif strcmpi(varargin{v},'Coordinate')
           Coordinate = varargin{v+1};
      elseif strcmpi(varargin{v},'Basics')
           Basics = varargin{v+1};
      elseif strcmpi(varargin{v},'Figure')
           Figure = varargin{v+1};
      end
   end
end

%% do the different jobs

disp('-----------------------------------')
fprintf(' 1st level QA - processing data %s \n',time_series);
disp('-----------------------------------')

if strcmpi(Basics,'on')
    [jobs.mean_image,jobs.std_image] = spmup_basics(time_series,'mean','std');
    V = spm_vol(jobs.mean_image);
else
    V = spm_vol(time_series);
end

% motion QA
if ~exist('Coordinate','var')
    Coordinate = V(1).dim /2;
end

[realign_files,jobs.FD,jobs.glo] = spmup_realign_qa(time_series,...
    'Motion Parameters', MotionParameters,...
    'Radius',Radius,...
    'Framewise Displacement', FramewiseDisplacement,...
    'Voltera', Voltera,...
    'Globals', Globals,...
    'Movie', Movie,...
    'Coordinate', Coordinate, ...
    'Figure',Figure);

if strcmp(FramewiseDisplacement,'on') || strcmp(Globals,'on')
    jobs.design = realign_files{1};
end

if strcmpi(Movie, 'on')
    jobs.movies = realign_files{2}; 
end

