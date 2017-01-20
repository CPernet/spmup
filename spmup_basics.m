function out = spmup_basics(P,operator)

% routine to compute a basic operation on data
%
% FORMAT out = spmup_basics(P,operator)
%
% INPUT P is a string of images to use
%       operator can be 'mean' 'std'
%
% OUTPUT if not specifiied write the result on drive
%        out the result from the specified operation
%
% Cyril Pernet Decembre 2016
% --------------------------------------------------------------------------
% Copyright (c) SPM Utility Plus toolbox


if nargin == 0
    [P,sts] = spm_select(Inf,'image','select the time series',[],pwd,'.*',1);
    if sts == 0 || isempty(P)
        disp('spmup_despike aborted'); return
    end
    operator = 'mean';
end
V = spm_vol(P);
Y = spm_read_vols(V);
dim = numel(size(Y));

if strcmp(operator,'mean')
    out = nanmean(Y,dim);
elseif strcmp(operator,'std')
    out = nanstd(Y,dim);
end

if nargout == 0
    [pathstr,name,ext]= fileparts(V(1).fname);
    V(1).fname = [pathstr filesep 'spmup_' operator '_' name ext];
    V(1).descrip = ['spmup ' operator];
    spm_write_vol(V(1),out);
end

