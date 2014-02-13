function spm_single_slices(varargin)

% this routine creates a series of single slices in a given orientation 
% with the fMRI map on top of it in a given color scheme
% each images are saved as jpg and also together as a .avi
%
% FORMAT spm_single_slices(anat,stat_map)
%        spm_single_slices(anat,stat_map,orientation,mycolor)
%
% INPUTS anat is the full name of the anatomical image
%        stat_map is the full name of the statistical image
%        mycolor is either a RGB vector or [];
%        orientation is 'axial', 'sagital' or 'coronal', 'all'
%                    if no input all three orientations are done
%
% OUTPUTS a series of jpg and a avi movie strored in a folder with the name
% of the stat map used
%
% Cyril Pernet 13-Feb-2014
% -------------------------
% Copyright (C) SPM-U+ 2014

%% input check
if nargin<2
    error('at least 2 imputs needed')
end
anat = varargin{1}; 
stat_map = varargin{2};
if ~ischar(anat) || ~ischar(stat_map)
    error('anat and stat_map must be names i.e. char')
end
    
[pathstr name ext] = fileparts(stat_map);
mycolor = hot(32);
orientation = 'all';
if nargin == 3
    mycolor = varargin{3};
    if ~isvector(mycolor)
        error('a [RGB] vector is expected for color')
    end
elseif nargin == 4
    mycolor = varargin{3};
    orientation = varargin{4};
end

%% quick data check
Vanat = spm_vol(anat);
Vstat_map = spm_vol(stat_map);
if sum(Vanat.dim == Vstat_map.dim) ~=3
    error('the anatomical and stat image have different dimensions, please resize')
end

%% load data and convert to 32 bits
% load the anatomical image and rescale for 32 bits
anat = spm_read_vols(Vanat);
scaled_anat = anat;
for z=1:size(anat,3)
    tmp = squeeze(anat(:,:,z));
    scaled_anat(:,:,z) =  round(31.*(tmp./ max(tmp(:))))+1;
end

% load the stat image and rescale for 32 bits
stat_map = spm_read_vols(Vstat_map);
scaled_stat_map = stat_map;
for z=1:size(stat_map,3)
    tmp = squeeze(stat_map(:,:,z));
    scaled_stat_map(:,:,z) =  round(31.*(tmp./ max(tmp(:))))+1;
end

%% colormap
% colormap in 64 bits
if isvector(mycolor)
    mycolor = repmat(mycolor,32,1);
end
mycolormap = [mycolor ; gray(32)];
   
%% figures
if strcmp(orientation,'all') || strcmp(orientation,'axial') || strcmp(orientation,'Axial')
    mkdir([name '_axial']); cd([name '_axial']);
    figure; colormap(mycolormap);
    vidObj = VideoWriter('axial.avi');
    open(vidObj); set(gca,'NextPlot','replacechildren')
    for a=1:Vanat.dim(3)
        image(squeeze(scaled_stat_map(:,:,a))'); 
        hold on; h=image(squeeze(scaled_anat(:,:,a))'+32); 
        set(h,'AlphaData',0.4); axis tight
        saveas(gcf,[name '_axial_' num2str(a) '.jpg'],'jpeg'); 
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    close(vidObj); close; 
    cd ..
end

if strcmp(orientation,'all') || strcmp(orientation,'sagital') || strcmp(orientation,'Sagital')
    mkdir([name '_sagital']); cd([name '_sagital']);
    figure; colormap(mycolormap);
    vidObj = VideoWriter('sagital.avi');
    open(vidObj); set(gca,'NextPlot','replacechildren')
    for a=1:Vanat.dim(1)
        image(squeeze(scaled_stat_map(a,:,:))'); 
        hold on; h=image(squeeze(scaled_anat(a,:,:))'+32); 
        set(h,'AlphaData',0.4); axis tight
        saveas(gcf,[name '_axial_' num2str(a) '.jpg'],'jpeg'); 
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    close(vidObj); close; 
    cd ..
end

if strcmp(orientation,'all') || strcmp(orientation,'coronal') || strcmp(orientation,'Coronal')
    mkdir([name '_coronal']); cd([name '_coronal']);
    figure; colormap(mycolormap);
    vidObj = VideoWriter('coronal.avi');
    open(vidObj); set(gca,'NextPlot','replacechildren')
    for a=1:Vanat.dim(3)
        image(squeeze(scaled_stat_map(:,a,:))'); 
        hold on; h=image(squeeze(scaled_anat(:,a,:))'+32); 
        set(h,'AlphaData',0.4); axis tight
        saveas(gcf,[name '_axial_' num2str(a) '.jpg'],'jpeg'); 
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    close(vidObj); close; 
    cd ..
end


    
    
