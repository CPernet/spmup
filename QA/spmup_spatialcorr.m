function [volume_outliers, slice_outliers] = spmup_spatialcorr(varargin)

% look at the time series per slice, and compute correlations
% this allows testing that images are alike and spots drop out easily
% also create the mean image, sd, and binary mask 
%
% FORMAT: [volume_outliers, slice_outliers] = spmup_spatialcorr
% FORMAT: [volume_outliers, slice_outliers] = spmup_spatialcorr(P,'figure',value)
%
% INPUT:  P are the fMRI data
%         'figure' is set as 'on' 'off' or 'save' (default)
%
% OUTPUT: volume_outliers indexes volumes with an average correlation between
%                         slices bigger or lower than the other volumes
%         slice_outliers indexes slices with an average correlation between
%                        volumes bigger or lower than the other slices
%
% Relies on SPM functions (see http://www.fil.ion.ucl.ac.uk/spm/)
% Cyril Pernet v4 - 25 Feb 2016, The University of Edinburgh
% ---------------------------------------------------------------

if exist('nansum','file') ~= 2
    error('you do not have stats toolbox to perform this operation, sorry')
end

%% check inputs

disp('running spmup_spatialcorr ...')
disp('-------------------------')

%% get data and mask
% memory mapped data
if nargin == 0
    [P,sts] = spm_select(Inf,'image' ,'Select your fMRI time series',{},pwd,'.*',Inf);
    if sts == 0
        return
    end
    V = spm_vol(P);
    % bypass orientation check allowing to look at raw data
    N = numel(V);
    Y = zeros([V(1).dim(1:3),N]);
    for i=1:N
        for p=1:V(1).dim(3)
            Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
        end
    end
else
    P = varargin{1};
    if ischar(P)
        V = spm_vol(P);
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
        
        for i=numel(V):-1:1
            for p=V(1).dim(3):-1:1
                Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
            end
        end
    else
        if numel(size(P)) == 4 % this is already data in
            Y = P; 
        else
            error('input data are not char nor 4D data matrix, please check inputs')
        end
    end
end
   
fig = 'save';
if nargin > 1 && strcmpi(varargin{2},'figure')
    fig = 'on';
    if nargin == 3
        fig = varargin{3};
    end
end

%% correlations between slices within each volume
% -----------------------------------------------

for j=size(Y,4):-1:1
    for i=(size(Y,3)-1):-1:1
        r(i,j) = corr2(Y(:,:,i,j),Y(:,:,i+1,j));
    end
end

if strcmpi(fig,'on') || strcmpi(fig,'save')
    figure('Name','Spatial QA')
    set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1); imagesc(r); ylabel('slices'); xlabel('volume nb')
    title('Correlations between slices within each volume')
end

r = nanmean(r,1); % average over slices
volume_outliers = spmup_comp_robust_outliers(r);
if sum(volume_outliers ~= 0)
    fprintf('volume %g is potiential outlier\n',find(volume_outliers));
else
    disp('no volume outliers detected');
end

if strcmpi(fig,'on') || strcmpi(fig,'save')
    subplot(2,2,3); plot(r,'LineWidth',2); ylabel('correlation value'); xlabel('volume nb')
    hold on; tmp = volume_outliers.*r; tmp(tmp==0)=NaN; plot(tmp,'ro','LineWidth',3);
    grid on; axis tight; hold on; title('Average correlations per volume'); clear tmp
end

%% correlations between volumes within each slice
% -----------------------------------------------

for j=(size(Y,4)-1):-1:1
    for i=size(Y,3):-1:1
        r(i,j) = corr2(Y(:,:,i,j),Y(:,:,i,j+1));
    end
end

if strcmpi(fig,'on') || strcmpi(fig,'save')
    subplot(2,2,2); imagesc(r)
    title('Correlations between volumes within each slice')
    ylabel('slices'); xlabel('volume nb')
end

r = nanmean(r,2); % average over volumes
slice_outliers = spmup_comp_robust_outliers(r);
if sum(slice_outliers ~= 0)
    fprintf('slice %g is potiential outlier\n',find(slice_outliers));
else
    disp('no slice outliers detected');
end

if strcmpi(fig,'on') || strcmpi(fig,'save')
    subplot(2,2,4); plot(r,'LIneWidth',2); ylabel('correlation value'); xlabel('slice nb')
    hold on; tmp = slice_outliers.*r; tmp(tmp==0)=NaN; plot(tmp,'ro','LineWidth',3);
    grid on; axis tight; hold on; title('Average correlations per slice'); clear tmp
end

if strcmpi(fig,'save')
    [filepath,filename] = fileparts(V(1).fname);
    try
        print (gcf,'-dpdf', '-bestfit', fullfile(filepath,[filename '_spatial_correlation.pdf']));
    catch
        print (gcf,'-dpdf', fullfile(filepath,[filename '_spatial_correlation.pdf']));
    end
    close(gcf)
end


