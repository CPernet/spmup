function [fmri_out,V4] = spmup_skip(varargin)

% routine that reads a 4D fmri data set and return a 4D set minus N initial volumes
%
% FORMAT fmri_out = spmup_skip
%        fmri_out = spmup_skip(fmri_in,n,options)
%
% INPUTS  fmri_in is the name of the fmri data set
%         n       is the nymber of initial volumes to skip
%         options 'bids_derivatives','on' (default) or 'off'  
%                 'rename','newname'        
%         - with bids_derivatives on there is no need of a new name, the same
%         name is used in /derivatives/sub-XXX/func assumimg your data are
%         already in bids (ie looking for a sub-xxx in the file name)
%         - rename is used to give a new name, next to fmri_in or in
%         derivatives depending if it's on or off
%
% OUTPUTS fmri_out is the name of the new reduced data set
%         V4 the spm_vol handle of the concatenated 3D files
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

if nargin == 0
    [P,sts] = spm_select(1,'image' ,'Select your fMRI time series',{},pwd,'.*',[1 Inf]);
    if sts == 0
        return
    end
    
    if strcmpi(P(end-1:end),',1')
        P = P(1:end-2);
    end
    
    V = spm_vol(P);
    N = numel(V);
    Y = zeros([V(1).dim(1:3),N]);
    for i=1:N
        for p=1:V(1).dim(3)
            Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
        end
    end
    
    skipped = inputdlg('how many initial volumes to skip?');
    if isempty(skipped)
        disp('spmup_skip aborded'); return
    else
       skipped = cell2mat(skipped);
    end
elseif nargin == 2
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
    end
    skipped = varargin{2};
end

%% reduce size
if ischar(skipped)
    skipped = eval(skipped);
end
Y(:,:,:,1:skipped) = [];

%% save
derivatives = 'on';
newname     = [];

if nargin > 2
    for n=3:nargin
        if contains(varargin{n},'derivatives')
            derivatives = varargin{v+1};
        elseif strcmpi(varargin{n},'rename')
            newname = varargin{n+1};
        end
    end
end

if strcmpi(derivatives,'on')
    index_name    = strfind(P,'sub-');
    index_name    = index_name(1); % just in case
    index_name(2) = strfind(P,[filesep 'func']);
    sub_name      = P(index_name(1):index_name(2)-1);
    root_name     = [P(1:index_name(1)-1) 'derivatives' filesep sub_name];
else
    [root_name,sub_name,ext] = fileparts(P);
end
    
if isempty(newname) 
    if strcmp(derivatives,'on')
        newname = fullfile(root_name, P(index_name(2):end));
    else
        newname = fullfile(root_name,['skiped_' sub_name ext]);
    end
else
    if ~strcmp(newname(end-3:end),'.nii')
        newname = [newname '.nii']; % if extension not provided in newname
    end
    
    if strcmp(derivatives,'on')
        newname = fullfile(root_name, ['func' filesep newname]);
    % else this is the input name = do nothing
    end
end

%% redo header
for i=(size(V,1)-skipped):-1:1
    W(i)       = V(i+skipped);
    W(i).fname = newname;
    W(i).n(1)  = i; % update which image in the series 
end

%% write to disk
if ~exist(fileparts(newname),'dir')
    mkdir(fileparts(newname))
end

for i=1:length(W)
    spm_write_vol(W(i),squeeze(Y(:,:,:,i)));
end
V4       = spm_file_merge(W,newname);
fmri_out = newname;

