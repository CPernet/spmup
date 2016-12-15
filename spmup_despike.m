function despiked = spmup_despike(varargin)
%
% SPM UP routine to 'despike' fMRI time-series in a similar way as AFNI does
% Note is requires the statistics toolbox (nansum, icdf are called)
%
% FORMAT spmup_despike
%        spmup_despike(P)
%        spmup_despike(P,M)
%        spmup_despike(P,M,flags)
%        spmup_despike(P,[],flags)
%
% INPUT if none the user is prompted
%       P the names of the fMRI images (time-series) or the 4D matrix of data
%       M the name of the mask or the 3D binary matrix
%       flags defines options to be used
%             - flags.auto_mask,'off' or 'on' if M is not provided, auto_mask is
%              'on' but if set to 'off' the user is prompted to select a mask
%             - flags.method is 'median' or any of the option of the 'smooth'
%                matlab function - in all cases the span is function of the
%                autocorrelation unless window is specified
%             - flags.window defines the number of consecutive images to use
%               to despike the data ; for instance flags.method = 'median'
%               and flags.window = 3 means that each data point is 1st
%               substituted by a moving 3 points median and the resulting fit
%               is used to determine outliers (see below)
%
% OUTPUT despiked is either the list of despiked images save onto disk or
%                 the despiked data, matching the input P
%        spmup_despike_log is saved onto disk where the data are
%        spmup_despike_log can be reviewed using spmup_review_despike_log
%        spmup_despike_log is structure with the fields:
%        - P the list of files (data) analyzed
%        - flags the flags used
%        - despiked_voxels the proportion of despiked voxels per volume
%        - class a 4d binary matrix indicating despiked voxels
%
%        the new dataset with the spikes removed is also written in the
%        data folder with the prefix 'despiked_' if no input or P is a
%        series of names given as input
%
% --------------------------------------------------------------------------
% The data are 1st smoothed either with a median filter (using the window
% parameter - see flag options) or using a smooth curve (see flag options)
% and then we look for outliers in the residuals and the data are
% interpolated. Although the smooting method is different, the detection of
% spikes and interpolation follows AFNI 3dDespike
% <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html>
% Finaly, data are written with the prefix despiked_  and a log file
% --------------------------------------------------------------------------
%
% Cyril Pernet Decembre 2016
% --------------------------------------------------------------------------
% Copyright (c) SPM Utility Plus toolbox

if exist('nansum','file') ~= 2
    error('you do not have stats toolbox to perform this operation, sorry')
end

if exist('smooth','file') ~= 2
    error('you need the curve fitting toolbox to perform this operation, sorry')
end

%% check inputs

% defaults
get_data = 1; % request data
get_mask = 0; % auto_mask
flags = struct('auto_mask','on','method','median');

% inputs
if nargin == 1
    get_data = 0;
    
elseif nargin == 2;
    get_data = 0;
    get_mask = 0;
    flags.auto_mask = 'off';
    
elseif nargin == 3
    get_data = 0;
    get_mask = 0;
    
    if isfield(varargin{3},'auto_mask')
        if strcmp(varargin{3}.auto_mask,'on') || strcmp(varargin{3}.auto_mask,'off')
            flags.auto_mask = varargin{3}.auto_mask;
        else
            error('flags.auto_mask must be ''on'' or ''off''');
        end
    end
    
    if isfield(varargin{3},'method')
        if strcmp(varargin{3}.method,'median')
            flags.method = varargin{3}.method;
        else
            if (exist('smooth','file') == 2) == 0 % ie doesn't exist
                error('flags.method must be ''median'' you don''t seem to have the curve fitting toolbox');
            elseif strcmp(varargin{3}.method,'moving') || strcmp(varargin{3}.method,'lowess') || ...
                    strcmp(varargin{3}.method,'loess') || strcmp(varargin{3}.method,'sgolay') || ...
                    strcmp(varargin{3}.method,'rlowess') || strcmp(varargin{3}.method,'rloess')
                flags.method = varargin{3}.method;
            else
                error('the method selected is not recognized ?? check flags.method');
            end
        end
    end
    
    if isfield(varargin{3},'window')
        flags.window = varargin{3}.window;
        if ~isempty(flags.window) && ~isnumeric(flags.window)
            error('the window parameter must be a real number')
        end
    end
end

disp('running spmup_despike ...')
disp('-------------------------')

%% get data and mask
% memory mapped data
if get_data == 1;
    [P,sts] = spm_select(Inf,'image','select the time series',[],pwd,'.*',1);
    V = spm_vol(P);
    % bypass orientation check
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
    else
        if numel(size(P)) == 4 % this is already data in
            Y = P; N = size(Y,4);
        else
            error('input data are not char nor 4D data matrix, please check inputs')
        end
    end
end


% memory mapped mask
if get_mask == 1;
    [M,sts] = spm_select(1,'image','select the mask',[],pwd,'.*',1);
    if sts ~=0
        Mask = spm_read_vols(spm_vol(M));
    else
        error('spm u+ stopped - mask selection interupted')
    end
elseif get_mask == 0 && strcmp(flags.auto_mask,'off')
    if ischar(varargin{2})
        try
            Mask = spm_read_vols(spm_vol(varargin{2}));
        catch varin_error
            error(varin_error.name)
        end
    else
        Mask = varargin{2};
        % just to make sure - binarize
        Mask(find(Mask)) = 1;
    end
else
    disp('generating a mask')
    Mask = spmup_auto_mask(V);
end

% figure('Name','Mask')
% colormap('gray')
% for z=1:size(Mask,3)
%     imagesc(flipud(Mask(:,:,z)'));
%     axis square; title(['Slice ' num2str(z)])
%     pause
% end

%% now do the despiking
% although the smooting method is different, the detection of spike and interpolation
% follows 3dDespike http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html

despiked_voxels = [];
if ~isfield(flags,'window') || isempty(flags.window)
    R = spmup_autocorrelation(V,Mask);
end

disp('smoothing data')
index = find(Mask);
YY = cell(1,length(index));
parfor i=1:length(index)
    [x,y,z]=ind2sub([size(Y,1),size(Y,2),size(Y,3)],index(1));
    data = squeeze(Y(x,y,z,:))';
    newdata = zeros(size(data));
    
    % define window size
    window = R(x,y,z); % need at least 3 points
    if window < 3; window = 3; end
    
    % median smoothing
    % ------------------
    if strcmp(flags.method,'median')
        
        % beginning
        for p=1:floor(window/2)
            newdata(p) = nanmedian([repmat(data(1),1,ceil(window/2)-p) data(1:p) data(p+1:p+floor(window/2))]);
        end
        % middle
        for p=(floor(window/2)+1):(N-floor(window/2))
            newdata(p) = nanmedian(data((p-floor(window/2)):(p+floor(window/2))));
        end
        % end
        last = 1;
        for p=(N-floor(window/2)+1):(N-1) % don't do last data point
            newdata(p) = nanmedian([data(p-ceil(window/2):p-1) repmat(data(p),1,floor(window/2)-last)]);
            last = last+1;
        end
        newdata(N) = data(N);
        
        % smooth function
        % --------------------------
    else
        newdata = smooth(data,window,flags.method);
    end
    
    % MAD of the residuals
    res = data-newdata;
    MAD = nanmedian(abs(res - repmat(nanmedian(res),1,N)));
    SIGMA = sqrt(pi/2)*MAD;
    s = res/SIGMA;
    
    %  * Values with s > c1 are replaced with a value that yields
    %     a modified s' = c1+(c2-c1)*tanh((s-c1)/(c2-c1)).
    %  * c1 is the threshold value of s for a 'spike' [default c1=2.5].
    %  * c2 is the upper range of the allowed deviation from the curve:
    %     s=[c1..infinity) is mapped to s'=[c1..c2)   [default c2=4].
    
    if SIGMA ~=0 % i.e. not res / 0
        out = find(s > 2.5);
        c1 = 2.5; c2=4; s2 = s;
        for p=1:length(out)
            s2(out(p)) = c1+(c2-c1)*tanh((s(out(p))-c1)/(c2-c1));
        end
        
        % reverse s2 to the real data
        YY{i} = (s2*SIGMA)+newdata;
    end
end

Despiked_QA = NaN(size(Y,1),size(Y,2),size(Y,3));
for i=1:index
    [x,y,z]=ind2sub([size(Y,1),size(Y,2),size(Y,3)],index(1));
    Despiked_QA(x,y,z) = sum(YY{i} ~= squeeze(Y(x,y,z,:))')/size(Y,4).*100;
    if ~isnan(Despiked_QA(x,y,z))
        Y(x,y,z,:) = YY{i};
    end
end
clear YY

    
%% write and return the data
if ischar(P)
    for v=1:size(Y,4)
        disp('writing data')
        V(v).descrip = 'spmup despiked';
        [pathstr,name,ext]= fileparts(V(v).fname);
        V(v).fname = [pathstr filesep 'despiked_' name ext];
        despiked{v} = V(v).fname;
        spm_write_vol(V(v),squeeze(YY(:,:,:,v)));
    end
    
    if size(P,1) == 1
        [pth,file,ext]=fileparts(P);
        V4 = spm_file_merge(V,[pth filesep 'despiked_' file ext]);
    end
else
    despiked = YY;
end
    
%% write the report
disp('saving spmup_despike_log')
if ischar(P)
    spmup_despike_log.P = P;
end
spmup_despike_log.flags = flags;
try
    if ~isfield(flags,'window') || isempty(spmup_despike_log.window)
        spmup_despike_log.window = R;
        V(1).descrip = 'spmup window size';
        [pathstr,name,ext]= fileparts(V(v).fname);
        V(1).fname = [pathstr filesep 'window' ext];
        spm_write_vol(V(1),W);
    else
        spmup_despike_log.window = window;
    end
end
spmup_despike_log.despiked_voxels = Despiked_QA;
save spmup_despike_log spmup_despike_log
disp('despiking done')
disp('--------------')
