function [meanimg,stdimg] = spmup_basics(P,varargin)

% routine to compute a basic mean and std on a trime series
%
% FORMAT out = spmup_basics(P,operators)
%
% INPUT P is a string of images to use
%       operators can be 'mean' 'std' (or both)
%
% OUTPUT name of files om the drive
%
% Cyril Pernet 
% --------------------------------------------------------------------------
% Copyright (c) SPM Utility Plus toolbox

meanimg = [];
stdimg  = [];

% check P
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
        error('input data are not char or cell, please check inputs')
    end
end
dim = numel(size(Y));

%% check operations to do
if nargin == 1
    operators{1} = 'mean';
    operators{2} = 'std';
else
    for v=nargin:-1:2
        operators{v-1} = varargin{v-1};
    end
end

%% compute and write image

[pathstr,name,ext]= fileparts(V(1).fname);
for v=1:length(operators)
    clear out
    if strcmpi(operators{v},'mean')
        out          = nanmean(Y,dim);
        if strcmpi(name(end-4:end),'_bold')
            V(1).fname   = [pathstr filesep name(1:end-5) '-mean_bold' ext];
        else
            V(1).fname   = [pathstr filesep name '_mean' ext];
        end
        V(1).descrip = ['spmup ' operators{v}];
        out          = spm_write_vol(V(1),out);
        meanimg      = out.fname;
    elseif strcmpi(operators{v},'std')
        try
            out = nanstd(Y,0,dim);
        catch
            warning('serious issue can''t get std! trying voxel-wise')    
            out = zeros(size(Y,1),size(Y,2),size(Y,3));
            for x=1:size(Y,1)
                for y=1:size(Y,2)
                    for z=1:size(Y,3)
                        try
                            out(x,y,z) = nanstd(squeeze(Y(x,y,z,:)),0);
                        catch
                            sprintf('serious issue voxel %g %g %g \n',x,y,z)
                        end
                    end
                end
            end
        end
        if strcmpi(name(end-4:end),'_bold')
            V(1).fname   = [pathstr filesep name(1:end-5) '-std_bold' ext];
        else
            V(1).fname   = [pathstr filesep 'std_' name ext];
        end
        V(1).descrip = ['spmup ' operators{v}];
        out          = spm_write_vol(V(1),out);
        stdimg       = out.fname;
    end
end


