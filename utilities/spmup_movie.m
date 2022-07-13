function out = spmup_movie(P,varargin)

% simple routine to create a movie of a time series
%
% FORMAT: spmup_plotmotion
%         spmup_plotmotion(data,'coordinates',vector,'filename',name)
%
% INPUTS: data can be the fmri name or actual data, if not providing, user is prompted
%         'coordinates' (optional) are the voxel positions to make movies around
%         'filename' (optional) is the name to use to save the file
%                    (default same as P, or pwd if P is already the data)
%         'showfig' (optional) set to 'on' of 'off' to display the movie
%
% OUTPUTs: avi file
%
% Cyril Pernet 
% ------------------------------
% Copyright (C) spmup team 

current = pwd;

%% check input
if nargin == 0
    [P,sts] = spm_select(Inf,'image' ,'Select your fMRI time series',{},pwd,'.*',Inf);
    if sts == 0
        return
    end
    V = spm_vol(P);
    [root,fname] = fileparts(V(1).fname);
    filename     = fullfile(root,fname);
    % bypass orientation check allowing to look at raw data
    N = numel(V);
    Y = zeros([V(1).dim(1:3),N]);
    for i=1:N
        for p=1:V(1).dim(3)
            Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
        end
    end
else
    if ischar(P)
        V = spm_vol(P);
        [root,fname] = fileparts(V(1).fname);
        filename     = fullfile(root,fname);
        N = numel(V);
        Y = zeros([V(1).dim(1:3),N]);
        for i=1:N
            for p=1:V(1).dim(3)
                Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
            end
        end
    elseif iscell(P)
        for v=size(P,1):-1:1
            V(v) =spm_vol(P{v});
        end
        [root,fname] = fileparts(V(1).fname);
        filename     = fullfile(root,fname);

        for i=numel(V):-1:1
            for p=V(1).dim(3):-1:1
                Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
            end
        end
    else
        if numel(size(P)) == 4 % this is already data in
            filename = fullfile(pwd,'time_series');
            Y = P;
        else
            error('input data are not char nor 4D data matrix, please check inputs')
        end
    end
end

n = size(Y,4);
if n == 1
    cd(current); 
    error('this is not a time series, only 1 image detected');
end

showfig = 'on';
coordinates = [];
if nargin >1
   for v=1:(nargin-1)
      if strcmpi(varargin{v},'filename')
          [root,fname] = fileparts(varargin{v+1});
          if isempty(root)
              root = pwd;
          end
          filename = fullfile(root,fname);
           
      elseif strcmpi(varargin{v},'coordinates')
           coordinates = varargin{v+1};
      elseif strcmpi(varargin{v},'showfig')
          showfig = varargin{v+1};
      end
   end
end

%% make movies
if isempty(coordinates)
    coordinates(1) = size(Y,1) / 2;
    coordinates(2) = size(Y,2) / 2;
    coordinates(3) = size(Y,3) / 2;
end
coordinates = round(coordinates);
x = squeeze(Y(coordinates(1),:,:,:));
y = squeeze(Y(:,coordinates(2),:,:));
z = squeeze(Y(:,:,coordinates(3),:));
clear Y

disp('making the time series movie ...')
vidObj = VideoWriter(filename);
open(vidObj);
figure('Name',filename,'Visible',showfig);
set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized', 'outerposition',[0 0 1 1]);
for i = 1:n
    clf
    subplot(1,3,1);
    imagesc(flipud(x(:,:,i)'));axis square
    title(['volume ' num2str(i)]);
    subplot(1,3,2);
    imagesc(flipud(y(:,:,i)'));axis square
    title(['volume ' num2str(i)]);
    subplot(1,3,3);
    imagesc(fliplr(z(:,:,i))');axis square
    title(['volume ' num2str(i)]);
    colormap('gray'); drawnow;
    CF = getframe(gcf);
    writeVideo(vidObj,CF);
end
close(vidObj); close
out = [filename '.avi'];
cd(current)
