function new_files = spmup_first_level_qa(varargin)

% routine calling spmup_realign_qa and spmup_normalize_qa
%
% INPUT new_files = spmup_first_level_qa (prompt user)
%       new_files = spmup_first_level_qa(T1_name,Images,flags)
%
%       T1_name is the full name of the T1 image to use as reference
%       Images is a cell array of the time series per session / subjects
%       flags define various options to be used in spmup_realign_qa and
%             spmup_normalize_qa
%       flags.motion_parameters = 'on' (default) or 'off'
%                   --> plot motion paramters and 1st derivatives
%                   --> generate new regressor.txt file with motion parameters and outliers detected using Carling's modified boxplot rule
%       flags.globals = 'on' (default) or 'off'
%                   --> compute and plot globals (see spm_globals)
%                   --> append the regressor.txt file with global outliers detected using Carling's modified boxplot rule
%       flags.volume_distance = 'off' (default) or 'on'
%                   --> compute the distance of each volume to the mean
%       flags.movie = 'off' (default) or 'on'
%       flags.AC = [46 64 37];
%                   --> generate three movies along the x, y, z planes defined based on the [0 0 0] point
%                   --> if flags.AC is empty, the user is asked to input this coordinate
%       flags.average = 'on' (default) or 'off'
%                   --> create the average of normalized images
%                   --> if flags.T1 is on, also plot the outlines of the
%                   difference with the normalized T1
%        flags.T1 = 'on' (default) or 'off'
%                   --> compute the difference between the T1 (assumed
%                       normalized) and the normalized images - plot the outlines
%
% Cyril Pernet February 2014
% v2 updated to work with the batch Sept 2015
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2015

% defaults
spm('defaults','FMRI')
global defaults
current = pwd;

%% T1 input
if nargin == 0
    T1 = spm_select(1,'image','select normalized T1 image');
else
    T1 = varargin{1};
end

%% Images input
folder = {};
if nargin <2
    go = 1;
    while go ~= 0
        [filenames,sts] = spm_select(Inf,'image','select normalized images');
        if sts == 0 || isempty(filenames)
            go = 0; break;
        else
            Normalized{go} = filenames;
            [folder{go},name,ext] = fileparts(Normalized{go}(1,:));
            fprintf('folder selected %s \n',folder{go});
            go = go+1;
        end
    end
else
    if iscell(varargin{2})
        for f=1:size(varargin{2},2)
            Normalized{f} = cell2mat(varargin{2}{f});
            [folder{f},name,ext] = fileparts(Normalized{f}(1,:));
        end
    else
        for f=1:size(varargin{2},1)
            Normalized{f} = varargin{2}(f,:);
            [folder{f},name,ext] = fileparts(Normalized{f}(1,:));
        end
    end
end

%% flag definition
flags = struct('motion_parameters','on','globals','on','volume_distance','off','movie','off', ...
    'AC', [], 'average','on', 'T1', 'on');
if isfield(varargin{3},'motion_parameters'); flags.motion_paramters = varargin{3}.motion_parameters; end
if isfield(varargin{3},'globals'); flags.globals = varargin{3}.globals; end
if isfield(varargin{3},'volume_distance'); flags.volume_distance = varargin{3}.volume_distance; end
if isfield(varargin{3},'movie'); flags.movie = varargin{3}.movie; end
if isfield(varargin{3},'T1'); flags.T1 = varargin{3}.T1; end
if isfield(varargin{3},'average'); flags.average = varargin{3}.average; end

%% iterate
for f=1:size(folder,2)
    cd(folder{f})
    fprintf('processing data in %s \n',folder{f});
    
    % motion test
    new_files{f} = spmup_realign_qa(Normalized{f},flags);
    
    % normalization test
    spmup_normalize_qa(flags,Normalized{f},T1);
end
cd(current)
disp('1st level QA done')
