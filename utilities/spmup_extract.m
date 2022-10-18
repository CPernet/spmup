function timeseries = spmup_extract(fmridata,graymatter,atlas,varargin)

% routine to extract timeseries data as per atlas
%
% FORMAT timeseries = spmup_extract(fmridata,gaymatter,atlas)
%        timeseries = spmup_extract(fmridata,gaymatter,atlas,LowPass,FreqSampling)
%
% INPUT fmridata file name of the fMRI time series volume m
%       gaymatter file name of a graymatter mask to use (must be same dimension as fMRI)
%       atlas filenale of an atlas to use (will be resampled if different size)
%       Optional LowPass: the bandpass boundary,FreqSampling = 1/TR 
%                if those parameters are guven, the matlab 'lowpass'
%                function is applied on extracted time series
%
% OUTPUT timeseries is a cell array of voxels for each ROI of the atlas intersecting
%                   with gray matter, dimension is X voxels * Y volumes
%
%        from there get correlations, means, etc, as needed, for instance:
%        avg = cellfun(@(x) mean(x), timeseries, 'UniformOutput', false)
%        con = corr(cell2mat(timeseries)');
%
% Cyril Pernet 
% --------------------------
%  Copyright (C) SPMUP Team 

if nargin == 0
    
    [fmridata, sts] = spm_select(1, 'image', 'Select fMRI data');
    if sts == 0
        disp('selection aborded');
        return
    end
    if size(fmridata, 1) == 1 && strcmp(fmridata(length(fmridata) - 1:end), ',1')
        fmridata = fmridata(1:length(fmridata) - 2); % in case picked 4D put left ,1
    end
    
    [graymatter, sts] = spm_select(1, 'image', 'Select grey matter tissue image');
    if sts == 0
        disp('selection aborted');
        return
    end
    
    [atlas, sts] = spm_select(1, 'image', 'Select atlas image image');
    if sts == 0
        disp('selection aborted');
        return
    end
end

% check fmridata is 4D
if iscell(fmridata)
    
    if strcmp(fmridata{1}(end - 1:end), ',1')
        fmridata{1} = fmridata{1}(1:end - 2);
    end
    
    if strcmp(fmridata{2}(end - 1:end), ',1')
        fmridata{2} = fmridata{2}(1:end - 2);
    end
    
    % get the data
    Vfmri = spm_vol(fmridata{1});
    
else % not a cell array
    if strcmp(fmridata(1, end - 1:end), ',1')
        fmridata = fmridata(1:end - 2);
    end
    Vfmri = spm_vol(fmridata(1, :));
end

disp   ('------------------------------')
fprintf(' running spmup_extract on \n%s\n',Vfmri(1).fname)
disp   ('------------------------------')

%% deal with the masks
if iscell(graymatter)
    graymatter = cell2mat(graymatter);
end

graymatter = spm_vol(graymatter);
if any(graymatter.dim ~= Vfmri(1).dim)
    error('Dimension issue between data and grey matter');
end

graymatter = spm_read_vols(spm_vol(graymatter));
if min(graymatter(graymatter~=0)) < 0.5
    warning('the grey matter mask has values as low as %g, you may consider thresholding it',min(graymatter(graymatter~=0)))
end

if iscell(atlas)
    atlas = cell2mat(atlas);
end

atlas = spm_vol(atlas);
if any(atlas.dim ~= Vfmri(1).dim)
    warning('Dimension issue between data and atlas - interpolating the atlas');
    
    clear matlabbatch;
    matlabbatch{1}.spm.util.bbox.image = {Vfmri(1).fname};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    bb = spm_jobman('run', matlabbatch);
    [~,atlas_filename,ext] = fileparts(atlas.fname);
    copyfile(atlas.fname,fullfile(fileparts(Vfmri(1).fname), [atlas_filename ext]));
    roi = fullfile(fileparts(Vfmri(1).fname), [atlas_filename ext]);
    atlas = spm_vol(cell2mat( ...
        spmup_resize(roi, bb{1}.bb, ...
        abs(diag(Vfmri(1).mat(1:3, 1:3)))')));
    delete(roi)
end
atlas = spm_read_vols(atlas);

labels = unique(uint8(atlas(uint8(atlas)~=0)));
warning('%g labels found, returning a cell array of the same size',length(labels));
timeseries = cell(length(labels),1);

for r = length(labels):-1:1
    [x, y, z]                        = ind2sub(size(atlas), intersect(find(atlas == r), find(graymatter)));
    tmp                              = spm_get_data(Vfmri, [x y z]')';
    tmp(find(sum(isnan(tmp), 2)), :) = []; %#ok<FNDSB> % removes rows of NaN;
    tmp(sum(tmp, 2) == 0, :)         = []; % remove 0
    timeseries{r}                    = tmp;
    clear tmp
end

if exist(fullfile(fileparts(Vfmri(1).fname), ['r' atlas_filename ext]),'file')
    delete(fullfile(fileparts(Vfmri(1).fname), ['r' atlas_filename ext]))
end

%% possibly low-pass filter data now
if nargin == 4
    lowpass = varargin{4};
    FreqSampling = varargin{5};
end

if exist('lowpass','var')
    for r = length(labels):-1:1
        timeseries{r} = lowpass(timeseries{r},LowPass,FreqSampling);
    end
end

