function [r_outliers,r_course] = spmup_volumecorr(varargin)

% compute the Pearson correlation between volumes of the time series 
% and return possible outliers among correlated volumes (this assumes of
% course that r is not biased by outlying voxels between volumes)
%
% FORMAT [r_course, outliers] = spmup_volumecorr
%        [r_course, outliers] = spmup_volumecorr(P)
%        [r_course, outliers] = spmup_volumecorr(P,'mask',M,'figure',value)
%
% INPUT P the fMRI data 
%       optionally input the 'mask' M 
%                            'figure' as 'on' 'off' or 'save' (default)
%
% OUTPUT r_course a vector of correlations between volumes
%        outliers 0/1 indicates which volumes are outliers
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

if exist('nansum','file') ~= 2
    error('you do not have stats toolbox to perform this operation, sorry')
end

%% check inputs

disp('running spmup_volumecorr ...')
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
    
M   = [];
fig = 'save';
if nargin > 1
    for in = 2:nargin
        if strcmpi(varargin{in},'figure')
            try
                fig = varargin{in+1};
            catch
                fig = 'on';
            end
            
        elseif strcmpi(varargin{in},'mask')
            M = varargin{in+1};
            if ischar(M)
                VM = spm_vol(M);
                Mask = spm_read_vols(VM);
            else
                if numel(size(M)) == 3 % this is already data in
                    Mask = m;
                else
                    error('mask data are not char nor 3D data matrix, please check inputs')
                end
            end
        end
    end
end

if isempty(M)
    disp('generating a mask')
    Mask = spmup_auto_mask(V);
end

% check dim
if any([size(Y,1) size(Y,2) size(Y,3)] ~= size(Mask))
    error('fmri data and mask have different dimentions %gx%gx%g vs %gx%gx%g',...
        size(Y,1),size(Y,2),size(Y,3),size(Mask,1),size(Mask,2),size(Mask,3))
end

%% compute the Pearson correlation  
data     = Y.*Mask; clear Y Mask; 
data     = reshape(data,[prod(V(1).dim),size(data,4)]);
data(sum(isnan(data),2) == size(data,2),:) = []; % remove empty voxels
r_course = corr(data); diag = 1:length(r_course)-1;
if strcmpi(fig,'on') || strcmpi(fig,'save')
    figure('Name','Volume correlation')
    set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    subplot(2,1,1); imagesc(r_course); title('Volume correlation matrix')
end  
r_course = [r_course(sub2ind(size(r_course),diag,(diag+1))) r_course(1,size(r_course,2))];

%% get the outliers from the r_course
r_outliers = spmup_comp_robust_outliers(r_course);
if strcmpi(fig,'on') || strcmpi(fig,'save')
    subplot(2,1,2); 
    plot(r_course,'LineWidth',2); grid on; title('Volume to volume correlation ');
    hold on; axis([1 size(r_course,2) min(r_course)-0.01*range(r_course) 1])
    tmp = r_course.*r_outliers; tmp(tmp==0) = NaN; plot(tmp,'ro','LineWidth',3)
    if strcmpi(fig,'save')
        [filepath,filename] = fileparts(V(1).fname);
        try
            print (gcf,'-dpdf', '-bestfit', fullfile(filepath,[filename '_volume_correlation.pdf']));
        catch
            print (gcf,'-dpdf', fullfile(filepath,[filename '_volume_correlation.pdf']));
        end
        close(gcf)
    end
end
fprintf('outlier volumes: %g\n',find(r_outliers));
fprintf('report: %s\n',fullfile(filepath,[filename '_volume_correlation.pdf']))
