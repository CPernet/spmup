function [M,MG] = spmup_timeseriesplot(fmridata, cn1, cn2, cn3, varargin)

% routine to produce plots a la Jonathan Power
%
% FORMAT spmup_timeseriesplot(fmridata, c1, c2, c3, options)
%        [M,MG] = spmup_timeseriesplot(fmridata, c1, c2, c3, ...
%                                                       'motion','on',...
%                                                       'nuisances','on',...
%                                                       'correlation','on', ...
%                                                       'figure', 'on')
%
% INPUT fmridata is a cell for the time series (see spm_select)
%                if a cell array with 2 elements, both are plotted
%       c1, c2, c3 are the tissue classes derived from the segmentation
%
%       several options are also available
%       'motion' 'on' (default) 'off' or a vector/matrix associated to the data
%               --> this represent the motion computed from the motion parameter
%              (see also spmup_FD.m)
%
%       'correlation' 'on' (default) 'off' or a vector/matrix associated to the data
%               --> for resting fMRI this is useful to check against motion
%               and nuisance that volumes don't get decorrelated
%               (see also spmup_volumecorr.m)
%
%       'nuisances' 'on' (default) 'off' or a matrix (n*2) of noise computed from c2 and c3
%               --> these typically reflect changes in globals, respiration and
%               cardiac cycles that are not well captured by motion correction
%               (see also spmup_nuisance.m)
%
%       'clean' is 'on' (default) or off to remove mean and linear trend
%       from the data 
%
%       'makefig' 'on' (default) 'off' to decide whether to plot the figure
%       or only return M which is the matrix containing the carpet plot.
%
%       Note that motion and correlation can be vectors or matrices, the
%       second columns or row is then assumed to mark outliers
%       All three cases can be computed automatically in this function or
%       using specific calls, if provided the matrices are plotted as such
%       and could therefore come from other computations - if a cell array
%       of two time series is provided, it assumes the same motion and nuisance
%       apply, taking data from the 1st set eg data before and after denoising
%
% OUTPUTS
%        M is a cell array of the top percentile voxels used to make the
%        voxplot M{1] is Grey matter, M{2} is white matter and M{3} csf
%        MG is the matrix of gray matter voxel intersecting c1 and Yeo's atlas
%        + 2 figures
%        - a voxplot showing all grey and white matter voxels
%          in time, associated to the traces in options
%        - corrrelations matrices as per Yeo 7 parcellation
%
% Example from subjects structure: 
%         spmup_timeseriesplot(subjects{1}.func{1},subjects{1}.tissues{1}{1},...
%              subjects{1}.tissues{1}{2},subjects{1}.tissues{1}{3},'figure', 'on')
%
% Reference: Power, J.D. (2016). A simple but useful way to assess fMRI
% scan qualities. NeuroImage
% <http://dx.doi.org/10.1016/j.neuroimage.2016.08.009>
%
% Cyril Pernet 
% --------------------------
%  Copyright (C) SPMUP Team 

%% check fMRI data in
if nargin == 0
    
    [fmridata, sts] = spm_select(1, 'image', 'Select fMRI data');
    if sts == 0
        disp('selection aborded');
        return
    end
    if size(fmridata, 1) == 1 && strcmp(fmridata(length(fmridata) - 1:end), ',1')
        fmridata = fmridata(1:length(fmridata) - 2); % in case picked 4D put left ,1
    end
    
    [cn1, sts] = spm_select(1, 'image', 'Select grey matter tissue image');
    if sts == 0
        disp('selection aborted');
        return
    end
    
    [cn2, sts] = spm_select(1, 'image', 'Select white matter tissue image');
    if sts == 0
        disp('selection aborted');
        return
    end
    
    [cn3, sts] = spm_select(1, 'image', 'Select csf tissue image');
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

filepath = fileparts(Vfmri(1).fname);
disp   ('------------------------------')
fprintf(' running spmup_timeseries on %s\n',Vfmri(1).fname)
disp   ('------------------------------')

%% deal with the masks
if iscell(cn1)
    cn1 = cell2mat(cn1);
end
greymatter = spm_vol(cn1);
if any(greymatter.dim ~= Vfmri(1).dim)
    disp('Dimension issue between data and grey matter');
end

if iscell(cn2)
    cn2 = cell2mat(cn2);
end
whitematter = spm_vol(cn2);
if any(whitematter.dim ~= Vfmri(1).dim)
    error('Dimension issue between data and white matter');
end

if iscell(cn3)
    cn3 = cell2mat(cn3);
end
csf = spm_vol(cn3);
if any(csf.dim ~= Vfmri(1).dim)
    error('Dimension issue between data and CSF');
end

c1 = spm_read_vols(spm_vol(greymatter));
c2 = spm_read_vols(spm_vol(whitematter));
c3 = spm_read_vols(spm_vol(csf));

% make each class is mutually exclusive by making voxel content is higher
% than the sum of the other + base threshold at 30%
threshold = 0.7;

c1  = c1 .* (c1 > (c2 + c3));
c2  = c2 .* (c2 > (c1 + c3));
c3  = c3 .* (c3 > (c1 + c2));
c1  = c1 .* (c1 > threshold);
c2  = c2 .* (c2 > threshold);
c3  = c3 .* (c3 > threshold);

%% options
motion      = 'on';
nuisances   = 'on';
correlation = 'on';
clean       = 'on'; 
makefig     = 'on';

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'motion')
        motion      = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'nuisances')
        nuisances   = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'correlation')
        correlation = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'figure')
        makefig = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'clean')
        clean = varargin{i + 1};
    end
end

% if no figure, no point computing extra stuff
if strcmpi(makefig, 'off')
    motion      = 'off';
    nuisances   = 'off';
    correlation = 'off';
end

figplot = 0;

% Motion
if strcmpi(motion, 'on') || ~exist('motion', 'var')
    rfile = dir([fileparts(Vfmri(1).fname) filesep '*.txt']);
    if isempty(rfile)
        [rfile, sts] = spm_select(1, 'txt', 'Select realignmnet parameters');
        if sts == 0
            disp('selection aborted');
            return
        end
    else
        if size(rfile ,1) > 1
            idx= find(arrayfun(@(x) strncmp(x.name, 'rp', 2), rfile));
            rfile = rfile(idx).name;
        else
            rfile = rfile.name;
        end
    end
    motion = spmup_FD([fileparts(Vfmri(1).fname) filesep rfile])';
    figplot = figplot + 1;
    
    % add outliers
    motion(2, :) = spmup_comp_robust_outliers(motion);
    
elseif isnumeric(motion)
    if size(motion, 1) == size(Vfmri, 1)
        motion = motion';
    end
    figplot = figplot + 1;
end

% Nuisance
if strcmpi(nuisances, 'on') || ~exist('nuisances', 'var')
    if iscell(fmridata)
        nuisances = spmup_nuisance(fmridata{1}, whitematter, csf);
    else
        nuisances = spmup_nuisance(fmridata, whitematter, csf);
    end
    nuisances = struct2array(nuisances)';
    figplot   = figplot + 1;
    
elseif isnumeric(nuisances)
    if size(nuisances, 1) == size(Vfmri, 1)
        nuisances = nuisances';
    end
    figplot = figplot + 1;
end

% Correlation
if strcmpi(correlation, 'on') || ~exist('correlation', 'var')
    if iscell(fmridata)
        [r_course(1, :), r_course(2, :)] = spmup_volumecorr(fmridata{1});
    else
        [r_course(1, :), r_course(2, :)] = spmup_volumecorr(fmridata);
    end
    figplot = figplot + 1;   
elseif isnumeric(correlation)
    if size(correlation, 1) == size(Vfmri, 1)
        r_course = correlation';
    end
    figplot = figplot + 1;
end

%% get matrix to plot from fMRI data
disp('Checking Yeo''s networks to organize GM voxels');
atlas_path     = fullfile(fileparts(mfilename('fullpath')), '..', 'external');
atlas_filename = 'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii';
roi            = fullfile(atlas_path, atlas_filename);

Vroi = spm_vol(roi);
if any(Vroi.dim ~= Vfmri(1).dim) % maybe we need to resize the ROI image
    
    % 1st check if there is a resampled one that match
    roi2 = fullfile(atlas_path, ['r' atlas_filename]);
    reslice_atlas = true;
    if exist(roi2, 'file')
        Vroi2 = spm_vol(roi2);
        if all(Vroi2.dim == Vfmri(1).dim)
            Vroi = Vroi2;
            reslice_atlas = false;
        end
        clear Vroi2;
    end
    
    if reslice_atlas
        % 2nd if not resample
        clear matlabbatch;
        matlabbatch{1}.spm.util.bbox.image = {Vfmri(1).fname};
        matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
        bb = spm_jobman('run', matlabbatch);
        copyfile(roi,fullfile(fileparts(Vfmri(1).fname), atlas_filename));
        roi = fullfile(fileparts(Vfmri(1).fname), atlas_filename);
        Vroi = spm_vol(cell2mat( ...
            spmup_resize(roi, bb{1}.bb, ...
            abs(diag(Vfmri(1).mat(1:3, 1:3)))')));
        delete(roi)
    end
end

% for each network, take gray matter voxels
ROIdata    = spm_read_vols(Vroi);
ROIdata    = round(ROIdata);
roi_values = unique(ROIdata);
roi_values(roi_values == 0) = [];

M = []; N = [];
for r = length(roi_values):-1:1
    [x, y, z]                        = ind2sub(size(ROIdata), intersect(find(ROIdata == r), find(c1)));
    tmp                              = spm_get_data(Vfmri, [x y z]')';
    tmp(find(sum(isnan(tmp), 2)), :) = []; %#ok<FNDSB> % removes rows of NaN;
    tmp(sum(tmp, 2) == 0, :)         = []; % remove 0
    M                                = [M; tmp];  
    N                                = [N; mean(tmp,1)];  % average of the ROI
    tmp                              = corr(tmp);  
    WM(r)                            = median(tmp(triu(true(size(tmp)),1))); % median upper triangle
    clear tmp
end
MG = M;

if exist(fullfile(fileparts(Vfmri(1).fname), ['r' atlas_filename]),'file')
    delete(fullfile(fileparts(Vfmri(1).fname), ['r' atlas_filename]))
end

line_index = size(M, 1);
NGM        = size(M, 1);

% get data for white matter and CSF
[x, y, z] = ind2sub(size(c3), [find(c2); find(c3)]);
tmp       = spm_get_data(Vfmri, [x y z]')';
mad_vox   = mad(tmp, 1, 2);
M         = [M; tmp(mad_vox >  prctile(mad_vox, round(NGM * 100 / size(tmp, 1) / 2)), :)];

% remove mean and linear trend
if strcmpi(clean,'on')
    X        = [linspace(0, 1, length(Vfmri))' ones(length(Vfmri), 1)];
    cleanM   = M - (X * (X \ M'))';
    clear M;
    M{1}     = cleanM;
else
    tmp      = M; clear M
    M{1}     = tmp; clear tmp
end

if iscell(fmridata) && size(fmridata, 2) == 2
    cleanM2 = spmup_timeseriesplot(fmridata{2}, ...
        cn1, ...
        cn2, ...
        cn3, ...
        'motion', 'off', ...
        'nuisance', 'off', ...
        'correlation', 'off', ...
        'makefig', 'off',...
        'clean', clean);
    M{2} = cleanM2;
    clear M2;
end

%% figure
if ~strcmpi(makefig, 'off')
    thin_line_width  = 2;
    thick_line_width = 2;
    figure('Name','Residuals Carpet Plot')
    if strcmpi(makefig,'on')
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    else
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1],'visible','off')
    end
    plotindex = 1;
    
    % top of the figure are nuisance time courses:
    % - head motion as framewise displacement
    % - the white mattrer and csf time courses (detrended)
    % - volume correlations
    
    if exist('motion', 'var')       
        subplot(figplot+4, 1, plotindex);
        if size(motion, 2) == 1
            motion = motion';
        end
        plot(motion(1, :), 'LineWidth', thick_line_width);
        plotindex = plotindex + 1;        
        hold on;
        if size(motion, 1) == 2
            tmp = motion(1, :) .* motion(2, :);
            tmp(tmp == 0) = NaN;
            plot(tmp, 'or', 'LineWidth', thin_line_width);
        end
        axis([1 length(motion) min(motion(1, :)) max(motion(1, :))]);
        ylabel('motion');
        title('Framewise Displacement');
        grid on;
    end
    
    if exist('nuisances', 'var')
        subplot(figplot + 4, 1, plotindex);
        plotindex = plotindex + 1;
        if size(nuisances,1) == 2
            nuisances = nuisances';
        end
        plot(nuisances, '--', 'LineWidth', thick_line_width);     
        axis([1 length(motion) min(nuisances(:)) max(nuisances(:))]);
        ylabel('mean');
        title('WM and CSF signals');
        grid on;       
    end
    
    if exist('r_course', 'var')
        subplot(figplot + 4, 1, plotindex);
        plotindex = plotindex + 1;
        plot(r_course(1, :), 'LineWidth', thick_line_width);
        hold on;
        if size(r_course, 1) == 2
            tmp = r_course(1, :) .* r_course(2, :);
            tmp(tmp == 0) = NaN;
            plot(tmp, 'or', 'LineWidth', thin_line_width);
        end
        axis([1 length(r_course) min(r_course(1, :)) max(r_course(1, :))]);
        ylabel('r');
        title('Volume Correlations');
        grid on;
    end
    
    % plot is organized as follow:
    %
    % - cortex (subdivided as Yeo et al networks)
    % - cerebellum
    % - nuclei
    %
    % // thick green line //
    %
    % white matter - shows top 30% most variables voxels
    % ventricules  - shows top 30% most variables voxels
    %
    
    if numel(M) == 1
        subplot(figplot + 4, 1, plotindex:figplot+4);        
        plot_carpet(M{1}, line_index, thin_line_width);
        ylabel('  CSF          White Matter           GM Networks');
    else
        subplot(figplot + 4, 1, plotindex:figplot+2);
        plot_carpet(M{1}, line_index, thin_line_width);
        ylabel('CSF WM GM');
        subplot(figplot + 4, 1, plotindex + 2:figplot + 4);
        plot_carpet(M{2}, line_index, thin_line_width);
        ylabel('CSF WM GM');
    end
    xlabel('Time (scans)');
    if strcmpi(makefig,'save')
        if exist(fullfile(filepath,'spmup_QC.ps'),'file')
            print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(filepath,'spmup_QC.ps'));
        else
            print (gcf,'-dpsc2', '-bestfit', fullfile(filepath,'spmup_QC.ps'));
        end
    end
    
    if ~strcmpi(makefig,'on')
        close(gcf);
    end
    
    % ----------------2nd figure for correlations only when explicit ----------------
    if strcmpi(makefig,'on')
        figure('Name','Correlation matrices')
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
        subplot(1,2,1); imagesc(corr(MG')); caxis([-1 1]); title('Correlations among all Gray Matter voxel time series')
        xlabel('Voxels'); ylabel('Voxels');
        subplot(1,2,2);  N = corr(N');
        N = N - eye(size(N)) + diag(WM); imagesc(N);
        caxis([-1 1]); title('Median corr within and corr between mean network time series')
        xlabel('Networks'); ylabel('Networks'); colorbar
    end
    
end

end

function plot_carpet(cleanM, line_index, thin_line_width)

imagesc(cleanM);
colormap('gray');
hold on;

% plot green line separating GM from the rest
plot( ...
    1:size(cleanM, 2), ...
    line_index .* ones(1, size(cleanM, 2)), ...
    'r', ...
    'LineWidth', thin_line_width);
end
