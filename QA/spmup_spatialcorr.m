function [volume_outliers, slice_outliers] = spmup_spatialcorr(varargin)

% look at the time series per slice, and compute correlations
% this allows testing that images are alike and spots drop out easily
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
% Cyril Pernet 
% --------------------------
%  Copyright (C) SPMUP Team 

if exist('nansum','file') ~= 2
    error('you do not have stats toolbox to perform this operation, sorry')
end

%% check inputs
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
            Y = P; V.fname = pwd;
        else
            error('input data are not char nor 4D data matrix, please check inputs')
        end
    end
end
   
% option
fig = 'save';
if nargin > 1 && strcmpi(varargin{2},'figure')
    fig = 'on';
    if nargin == 3
        fig = varargin{3};
    end
end

% filename
if ischar(P)
    [filepath,filename] = fileparts(V(1).fname); 
    disp   ('-------------------------------')
    fprintf('running spmup_spatialcorr on %s\n',filename)
    disp   ('-------------------------------')
else
    disp   ('-------------------------')
    disp   ('running spmup_spatialcorr')
    disp   ('-------------------------')
end
filename(strfind(filename,'_')) = ' ';


%% correlations between slices for each volume
% -----------------------------------------------
n = size(Y,3);
for i = size(Y,4):-1:1
    S      = squeeze(Y(:,:,:,i));                % take one volume  
    vols   = reshape(S,[size(Y,1)*size(Y,2),n]); % reshape each slice as vectors
    all_r  = corr(vols);                         % all correlations between slices
    idx    = 2:n+1:numel(all_r);                 % take only one above diag (12,23,34,..) 
    r(i,:) = all_r(idx);
end

if strcmpi(fig,'on') || strcmpi(fig,'save')
    figure('Name','Spatial QA')
    if strcmpi(fig,'on')
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    else
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1],'visible','off')
    end
    subplot(2,2,1); imagesc(r'); ylabel('slices'); xlabel('volume nb')
    title(sprintf('Correlations between slices \n within each volume'))
end

r = nanmean(r,2); % average over slices
volume_outliers = spmup_comp_robust_outliers(r);
if sum(volume_outliers ~= 0)
    fprintf('volume %g is potiential outlier\n',find(volume_outliers));
else
    disp('no volume outliers detected');
end

if strcmpi(fig,'on') || strcmpi(fig,'save')
    subplot(2,2,2); plot(r,'LineWidth',2); ylabel('correlation value'); xlabel('volume nb')
    hold on; tmp = volume_outliers.*r; tmp(tmp==0)=NaN; plot(tmp,'ro','LineWidth',3);
    grid on; axis tight; hold on; title(sprintf('Correlations averaged across slices \n %s',filename)); 
    clear tmp
end
clear S vols all_r idx r

%% correlations between volumes slice by slice
% -----------------------------------------------

n = size(Y,4);
for i = size(Y,3):-1:1
    S      = squeeze(Y(:,:,i,:));                % take one slice all volumes
    vols   = reshape(S,[size(Y,1)*size(Y,2),n]); % reshape slices as vectors
    all_r  = corr(vols);                         % all correlations between volumes
    idx    = 2:n+1:numel(all_r);                 % take only one above diag (12,23,34,..) 
    r(i,:) = all_r(idx);
end

if strcmpi(fig,'on') || strcmpi(fig,'save')
    subplot(2,2,3); imagesc(r')
    title(sprintf('Correlations between volumes \n for each slice'))
    xlabel('slices'); ylabel('volume nb')
end

r = nanmean(r,2); % average over volumes
tmp = r(~isnan(r)); % if 0 in all slices, r is nan
tmp = spmup_comp_robust_outliers(tmp);
slice_outliers = zeros(length(r),1);
slice_outliers(~isnan(r)) = tmp;
if sum(slice_outliers ~= 0)
    fprintf('slice %g is potiential outlier\n',find(slice_outliers));
else
    disp('no slice outliers detected');
end

if strcmpi(fig,'on') || strcmpi(fig,'save')
    subplot(2,2,4); plot(r,'LineWidth',2); ylabel('correlation value'); xlabel('slice nb')
    hold on; tmp = slice_outliers.*r; tmp(tmp==0)=NaN; plot(tmp,'ro','LineWidth',3);
    grid on; axis tight; hold on; title(sprintf('Correlations averaged across volumes \n %s',filename)); 
    clear tmp
end

if strcmpi(fig,'save')
    if exist(fullfile(filepath,'spmup_QC.ps'),'file')
        print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(filepath,'spmup_QC.ps'));
    else
        print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(filepath,'spmup_QC.ps'));
    end
    close(gcf)
end


