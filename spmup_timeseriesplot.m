function M = spmup_timeseriesplot(fmridata,c1,c2,c3,varargin)

% routine to produce plots a la Jonathan Power
%
% FORMAT spmup_timeseriesplot(P,c1, c2, c3, options)
%        M = spmup_timeseriesplot(P,c1, c2, c3, 'motion','on','nuisances','on','correlation','on')
%
% INPUT fmridata is a cell array for the time series (see spm_select)
%       c1, c2, c3 are the tissue classes derived from the segmentation
%       
%       several options are also available
%       'motion' 'on' (default) 'off' or a vector/matrix associated to the data
%               --> this represent the motion computed from the motion parameter
%              (see also spmup_FD.m)
%       'correlation' 'on' (default) 'off' or a vector/matrix associated to the data
%               --> for resting fMRI this is useful to check against motion
%               and nuisance that volumes don't get decorrelated 
%               (see also spmup_volumecorr.m)
%       'nuisances' 'on' (default) 'off' or a matrix (n*2) of noise computed from c2 and c3
%               --> these typically reflect changes in globals, respiration and
%               cardiac cycles that are not well captured by motion correction
%               (see also spmup_nuisance.m)
%
%       Note that motion and correlation can be vectors or matrices, the
%       second columns or row is then assumed to mark outliers
%       All three cases can be computed automatically in this function or
%       using specific calls, if provided the matrices are plotted as such
%       and could therefore come from other computations
%
% OUTPUT a figure (voxplot) showing all grey and white matter voxels
%        in time, associated to the traces in options
%        M is the matrix of voxel by time
%
% Reference: Power, J.D. (2016). A simple but useful way to assess fMRI
% scan qualities. NeuroImage
% <http://dx.doi.org/10.1016/j.neuroimage.2016.08.009>
%       
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

%% check fMRI data in
if nargin == 0
    [fmridata,sts]= spm_select(Inf,'image','Select fMRI data');
    if sts == 0; disp('selection aborded'); return; end
    if size(fmridata,1) == 1 && strcmp(fmridata(length(fmridata)-1:end),',1')
        fmridata = fmridata(1:length(fmridata)-2); % in case picked 4D put left ,1
    end
    [c1,sts]= spm_select(1,'image','Select grey matter tissue image');
    if sts == 0; disp('selection aborded'); return; end
    [c2,sts]= spm_select(1,'image','Select white matter tissue image');
    if sts == 0; disp('selection aborded'); return; end
    [c3,sts]= spm_select(1,'image','Select csf tissue image');
    if sts == 0; disp('selection aborded'); return; end
end

% check fmridata is 4D
if iscell(fmridata)
    for v=1:size(fmridata,1)
        Vfmri(v) =spm_vol(fmridata{v});
    end
else
    Vfmri = spm_vol(fmridata);
end

if sum(size(Vfmri)) == 2
    error('fMRI data must be time series')
end

%% deal with the masks
if iscell(c1); c1=cell2mat(c1); end
greymatter = spm_vol(c1); 
if any(greymatter.dim ~= Vfmri(1).dim)
    disp('Dimension issue between data and grey matter - reslicing ... ')
end

if iscell(c2); c2=cell2mat(c2); end
whitematter = spm_vol(c2); 
if any(whitematter.dim ~= Vfmri(1).dim)
    error('Dimension issue between data and white matter')
end

if iscell(c3); c3=cell2mat(c3); end
csf = spm_vol(c3); 
if any(csf.dim ~= Vfmri(1).dim)
    error('Dimension issue between data and CSF')
end

c1 = spm_read_vols(spm_vol(greymatter));
c2 = spm_read_vols(spm_vol(whitematter));
c3 = spm_read_vols(spm_vol(csf));

% make each class is mutually exclusive by making voxel content is higher
% than the sum of the other + base threshold at 20%
c1  = c1 .* (c1>(c2+c3));  
c2  = c2 .* (c2>(c1+c3));  
c3  = c3 .* (c3>(c1+c2));  
c1  = c1 .* (c1>0.7);
c2  = c2 .* (c2>0.7);
c3  = c3 .* (c3>0.7);


motion = 'on';
nuisances = 'on';
correlation = 'on';

% options
for i=1:length(varargin)
   if strcmpi(varargin{i},'motion') 
       motion = varargin{i+1};
   elseif strcmpi(varargin{i},'nuisances') 
       nuisances = varargin{i+1};
   elseif strcmpi(varargin{i},'correlation') 
       correlation = varargin{i+1};
   end 
end

figplot = 0;
% Motion
if strcmpi(motion,'on') || ~exist('motion','var')
       rfile = dir([fileparts(Vfmri(1).fname) filesep '*.txt']);
    if isempty(rfile)
        [rfile,sts]= spm_select(1,'txt','Select realignmnet parameters');
        if sts == 0; disp('selection aborded'); return; end
    else
        rfile = rfile.name; 
    end
    motion = spmup_FD([fileparts(Vfmri(1).fname) filesep rfile])'; 
    figplot = figplot+1;
    
    % add outliers
    motion(2,:) = spmup_comp_robust_outliers(motion);

elseif isnumeric(motion)
    motion = motion; 
    if size(motion,1) == size(Vfmri,1)
        motion = motion';
    end
    figplot = figplot+1;
end

% Nuisance
if strcmpi(nuisances,'on') || ~exist('nuisances','var')
    nuisances = spmup_nuisance(fmridata,whitematter,csf); 
    nuisances = struct2array(nuisances)';
    figplot = figplot+1;
    
elseif isnumeric(nuisances) 
    if size(nuisances,1) == size(Vfmri,1)
        nuisances = nuisances';
    end
    figplot = figplot+1;
end

% Correlation
if strcmpi(correlation,'on') || ~exist('correlation','var')
    [r_course(1,:), r_course(2,:)] = spmup_volumecorr(fmridata);
    figplot = figplot+1;
    
elseif isnumeric(correlation) 
    if size(correlation,1) == size(Vfmri,1)
        r_course = correlation';
    end
    figplot = figplot+1;
end

% get matrix to plot from fMRI data
M = [];
roi = [fileparts(which('spmup_timeseriesplot.m')) filesep 'external' filesep 'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'];
Vroi = spm_vol(roi);
if any(Vroi.dim ~= Vfmri(1).dim) % maybe we need to resize the ROI image
    clear matlabbatch
    matlabbatch{1}.spm.util.bbox.image = {Vfmri(1).fname};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    bb = spm_jobman('run',matlabbatch);
    Vroi = spm_vol(cell2mat(spmup_resize(roi,bb{1}.bb,abs(diag(Vfmri(1).mat(1:3,1:3)))')));
end

% for each network, take gray matter voxels
ROIdata = spm_read_vols(Vroi);
roi_values = unique(ROIdata);
roi_values(roi_values==0) = [];
for r = 1:length(roi_values)
    [x,y,z]=ind2sub(size(ROIdata),intersect(find(ROIdata==r),find(c1)));
    tmp = spm_get_data(Vfmri,[x y z]')';
    tmp(find(sum(isnan(tmp),2)),:) = []; % removes rows of NaN;
    M = [M ; tmp];
end
line_index = size(M,1);

% get data for white matter and CSF
[x,y,z]=ind2sub(size(c3),[find(c2);find(c3)]);
M = [M ; spm_get_data(Vfmri,[x y z]')'];

% remove mean and linear trend
X = [linspace(0,1,length(Vfmri))' ones(length(Vfmri),1)]; 
cleanM = M - (X*(X\M'))'; 



%% figure

figure('Name','voxplot')
set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
plotindex = 1;
% top of the figure are nuisance time courses: head motion as framewise
% displacement and the white mattrer and csf time courses (detrended)
if exist('motion','var')
    subplot(figplot+4,1,plotindex); plotindex = plotindex+1;
    plot(motion(1,:),'LineWidth',3); hold on;
    if size(motion,1) == 2
        plot(motion(1,:).*motion(2,:),'or','LineWidth',2);
    end
    axis([1 length(motion) 1/10*range(motion(1,:)) max(motion(1,:))])
    ylabel('motion'); title('Framewise Displacement'); grid on; 
end

if exist('nuisances','var')
    subplot(figplot+4,1,plotindex); plotindex = plotindex+1;
    plot(nuisances','LineWidth',3); hold on;
    axis([1 length(motion) min(nuisances(:)) max(nuisances(:))])
    ylabel('mean'); title('WM and CSF signals'); grid on; 
end

% displacement and the white mattrer and csf time courses (detrended)
if exist('r_course','var')
    subplot(figplot+4,1,plotindex); plotindex = plotindex+1;
    plot(r_course(1,:),'LineWidth',3); hold on;
    if size(r_course,1) == 2
        plot(r_course(1,:).*r_course(2,:),'or','LineWidth',2);
    end
    axis([1 length(r_course) min(r_course(1,:)) max(r_course(1,:))])
    ylabel('r'); title('Volume Correlations'); grid on; 
end

% plot is organized as follow: cortex (subdivided as Yeo et al networks), cerebellum, nuclei 
% // thick line // white matter, ventricules
subplot(figplot+4,1,[plotindex:figplot+4]);
imagesc(cleanM); 
colormap('gray'); 
hold on
plot([1:size(M,2)],line_index.*ones(1,size(M,2)),'g','LineWidth',2)
xlabel('Time (scans)'); ylabel('  CSF          White Matter           GM Networks');


saveas(gcf, fullfile(fileparts(fmridata), 'voxplot.fig'),'fig'); 
try
    print (gcf,'-dpsc2', '-bestfit', fullfile(fileparts(fmridata), 'voxplot.ps'));
catch
    print (gcf,'-dpsc2', fullfile(fileparts(fmridata), 'voxplot.ps'));
end
close(gcf)

