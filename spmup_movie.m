function spmup_movie(data)

% simple routine to create a movie of a time series
%
% FORMAT: spmup_plotmotion(data)
%
% INPUT: the data to make the movie from
%        in not providing, user is prompted
%
% OUTPUT: a figure (both .fig and .jpg)
%
% Note: if you input raw data, it is likely that spm_read_vols complains
% 'Error using spm_check_orientations (line 65)'
% 'The orientations etc must be identical for this procedure.'
% 'Error in spm_read_vols (line 26)'
% 'spm_check_orientations(V);'
% Make sense - simply comment temporarily line 26 in spm_read_vols
% 
% Cyril Pernet - 02 Dec 2016
% ------------------------------
% Copyright (C) spmup team 2016

%% get the data
current = pwd;

if nargin == 1
    V = spm_vol(data);
else
    [P,sts]=spm_select(Inf,image,'select time series');
    V = spm_vol(P);
end

if size(V,1) == 1
    cd(current); error('this is not a time series, only 1 image detected');
end

%% make movies
Images = spm_read_vols(V); n = size(Images,4);
AC(1) = round(size(Images,1) / 2);
AC(2) = round(size(Images,2) / 2);
AC(3) = round(size(Images,3) / 2);
X = squeeze(Images(AC(1),:,:,:));
Y = squeeze(Images(:,AC(2),:,:));
Z = squeeze(Images(:,:,AC(3),:));
clear Images
    
plane = questdlg('Which orientation to plot data?', 'Plane question', ...
    'Axial', 'Coronal', 'Sagital', 'axial');

figure; colormap('gray')

% make video object
if strcmpi(plane,'sagital')
    vidObj = VideoWriter('sagital.avi');
elseif strcmpi(plane,'coronal')
    vidObj = VideoWriter('coronal.avi');
elseif strcmpi(plane,'axial')
    vidObj = VideoWriter('axial.avi');
end
open(vidObj);
set(gca,'NextPlot','replacechildren')

% make pictures and get frames
for i=1:n
    if strcmpi(plane,'sagital')
        imagesc(X(:,:,i)');axis tight
    elseif strcmpi(plane,'coronal')
        imagesc(Y(:,:,i)');axis tight
    elseif strcmpi(plane,'axial')
        imagesc(Z(:,:,i)');axis tight
    end
    title(['volume ' num2str(i)]);
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
end
close(vidObj);close

cd(current)
