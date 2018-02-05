function [despiked,filtered] = spmup_despike(varargin)
%
% SPM UP routine to 'despike' fMRI time-series.
% Note is requires the statistics toolbox (nansum, icdf are called)
%
% FORMAT spmup_despike
%        spmup_despike(P) 
%        spmup_despike(P,M) 
%        spmup_despike(P,M,flags)
%        spmup_despike(P,[],flags)
%
% INPUT if none the user is prompted
%       P the names of the 3D or 4D fMRI images (time-series) 
%       M the name of the mask or the 3D binary matrix 
%       flags defines options to be used
%             - flags.auto_mask,'off' or 'on' if M is not provided,
%             auto_mask is 'on' but if set to 'off' the user is prompted 
%             to select a mask
%             - flags.method is 'median' or any of the option of the
%             'smooth' matlab function - in all cases the span is function 
%             of the autocorrelation unless window is specified
%             - flags.window defines the number of consecutive images to
%             use to despike the data ; for instance flags.method = 'median'
%             and flags.window = 3 means that each data point is 1st
%             substituted by a moving 3 points median and the resulting
%             fit is used to determine outliers (see below)
%             - flags.skip defines the number of initial images to skip
%               often despiking is the 1st step in preprocessing and
%               therefore some initial scans needs to be discarded because
%               the series hasn't reach staty-state (unless prescan
%               normalizalise)
%
% OUTPUT despiked is either the list of despiked images save onto disk or
%                 the despiked data, matching the input P
%        filtered is a data median smoothed from which one despikes
%
%        spmup_despike_log is saved onto disk where the data are
%        spmup_despike_log can be reviewed using spmup_review_despike_log
%        spmup_despike_log is structure with the fields: 
%        - P the file name of the data despiked
%        - flags the option chosen
%        - window can be as flags or the full autocorrelation matrix
%        - RMS the root mean square difference netween the original and despiked data
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
%
% The prefered method is 'median' and window = []; this implies that the
% window size if based on the autocorrelation rather than fixed and the
% computation is fast because performed per window size - by contrast
% smoother option rely on estimates per voxel and this is much slower
% --------------------------------------------------------------------------
%
% Cyril Pernet Decembre 2016
% --------------------------------------------------------------------------
% Copyright (c) SPM Utility Plus toolbox

if exist('nansum','file') ~= 2
    error('you do not have stats toolbox to perform this operation, sorry')
end

%% check inputs

% defaults
get_data = 1; % request data
get_mask = 0; % auto_mask
flags = struct('auto_mask','on','method','median');

% inputs
if nargin == 1
    get_data = 0;
    
elseif nargin == 2
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
    
    if isfield(varargin{3},'skip')
        flags.skip = varargin{3}.skip;
        if ~isempty(flags.skip) && ~isnumeric(flags.skip)
            error('the skip parameter must be a real number')
        end
    else
        flags.skip = 0;
    end
end

disp('running spmup_despike ...')
disp('-------------------------')

%% get data and mask
% memory mapped data
if get_data == 1
    [P,sts] = spm_select(Inf,'image','select the time series',[],pwd,'.*',1);
    if sts == 0 || isempty(P)
        disp('spmup_despike aborted'); return
    end
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
        V(1:flags.skip) = []; % remove initial volumes
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
if get_mask == 1
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
    Mask = spmup_auto_mask(V());
end


%% now do the despiking
% although the smooting method is different, the detection of spike and
% interpolation follows 3dDespike
% http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html

if ~isfield(flags,'window') || isempty(flags.window)
    R = spmup_autocorrelation(V,Mask); 
    % figure; for z=1:size(R,3);imagesc(flipud(squeeze(R(:,:,z))')); pause; end
else
    R = ones(size(Y,1),size(Y,2),size(Y,3)).*flags.window;
end

disp('smoothing data')
if strcmp(flags.method,'median')
    % doesn't matter how fast is the autocorr coef - min is 3 TR to smooth
    index = unique(R(:)); index(index<=2) = 3;      % minimum window size is 3
    index(isnan(index))=[]; index = unique(index);  % the different window sizes to use
    if nargout == 2; filtered = NaN(size(Y)); end   % if we want to save filtered data
    YY = NaN(size(Y));                              % that's the despiked data
   
    for i=1:length(index)
        if index(i) == 3
            [x,y,z]=ind2sub([size(Y,1),size(Y,2),size(Y,3)],find(R<=index(i))); % which voxels to do
        else
            [x,y,z]=ind2sub([size(Y,1),size(Y,2),size(Y,3)],find(R==index(i)));
        end
        
        if ischar(P)
            data = spm_get_data(V,[x y z]');
            data = data'; 
        else
            data = NaN(length(x),N);
            for coord = 1:length(x)
                data(coord,:) = squeeze(Y(x(coord),y(coord),z(coord),:));
            end
        end
        newdata = zeros(size(data)); % this is the filtered data

        % beginning
        for p=1:floor(index(i)/2)
            newdata(:,p) = nanmedian([repmat(data(:,1),1,ceil(index(i)/2)-p) data(:,1:p) ...
                data(:,p+1:p+floor(index(i)/2))],2); % repeate half previous points, current point, next set of points
        end
        % middle
        for p=(floor(index(i)/2)+1):(N-floor(index(i)/2))
            newdata(:,p) = nanmedian(data(:,(p-floor(index(i)/2)):(p+floor(index(i)/2))),2);
        end
        % end
        last = 1;
        for p=(N-floor(index(i)/2)+1):(N-1) % don't do last data point
            newdata(:,p) = nanmedian([data(:,p-ceil(index(i)/2):p-1) repmat(data(:,p),1,floor(index(i)/2)-last)],2);
            last = last+1;
        end
        newdata(:,N) = data(:,N);
        
        % MAD of the residuals
        res = data-newdata;
        MAD = nanmedian(abs(res - repmat(nanmedian(res,2),1,N)),2);
        SIGMA = sqrt(pi/2).*MAD;
        s = res./repmat(SIGMA,1,N);
        
        %  * Values with s > c1 are replaced with a value that yields
        %     a modified s' = c1+(c2-c1)*tanh((s-c1)/(c2-c1)).
        %  * c1 is the threshold value of s for a 'spike' [default c1=2.5].
        %  * c2 is the upper range of the allowed deviation from the curve:
        %     s=[c1..infinity) is mapped to s'=[c1..c2)   [default c2=4].
        
        c1 = 2.5; c2=4; s2 = s;

        out = find(s > c1); % positive outliers
        for p=1:length(out)
            s2(out(p)) = c1+(c2-c1)*tanh((s(out(p))-c1)/(c2-c1));
        end
        
        out = find(s < -c1); % negative outliers
        for p=1:length(out)
            s2(out(p)) = -c1+(c2-c1)*tanh((s(out(p))+c1)/(c2-c1));
        end
        
        % reverse s2 to the real data ie multiply by var and add back the
        % smooth data (quick check: it works for test = (s.*repmat(SIGMA,1,N))+newdata;
        % sum(single(test(:)) == data(:)) is numel(data) to the rounding precision
        s2 = (s2.*repmat(SIGMA,1,N)) + newdata;
        
        % we can use indices here and remove that coord loop ?
        for coord = 1:size(x,1)
            if nargout == 2
                filtered(x(coord),y(coord),z(coord),:) = newdata(coord,:);
            end
            
            YY(x(coord),y(coord),z(coord),:) = s2(coord,:);
            % when SIGMA = 0 we have NaN -> back to original
            if isnan(s2(coord,1))
                YY(x(coord),y(coord),z(coord),:) = Y(x(coord),y(coord),z(coord),:);
            end
        end
    end
    
else  % smooth function
    disp('using matlab smoother ... this iterate per voxel and takes lots of time')
    
    if isempty(Mask)
        index = find(squeeze(Y(:,:,:,end)));
    else
        index = find(Mask);
    end
    YY = cell(1,length(index));
    ZZ = cell(1,length(index));
    
    parpool('local',feature('numCores')-1); % use all available cores -1
    parfor i=1:length(index)
        [x,y,z]=ind2sub([size(Y,1),size(Y,2),size(Y,3)],index(i));
        data = squeeze(Y(x,y,z,:))';
        newdata = zeros(size(data));
        
        % define window size
        window = R(x,y,z); % need at least 3 points
        if window < 3; window = 3; end
        if strcmp(flags.method,'moving')
            newdata = smooth(data,window)';
        else
            newdata = smooth(data,window,flags.method)';
        end
        ZZ{i} = newdata;
        
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
        else
            YY{i} = data;
        end
    end
    parpool close
end


% quick cleanup
if ~isempty(Mask)
    for v=1:size(Y,4)
        Y(:,:,:,v)  = Y(:,:,:,v).*Mask;
        YY(:,:,:,v) = YY(:,:,:,v).*Mask;
    end
end

%% write and return the data
if ischar(P)
    disp('writing data')
    if size(P,1) == size(YY,4) % 3D
        for v=1:size(Y,4)
            V(v).descrip = 'spmup despiked';
            [pathstr,name,ext]= fileparts(V(v).fname);
            V(v).fname = [pathstr filesep 'despiked_' name ext];
            despiked{v} = V(v).fname;
            spm_write_vol(V(v),squeeze(YY(:,:,:,v)));
        end
    else % 4D
        [pathstr,name,ext]= fileparts(V(1).fname);
        fname = [pathstr filesep 'despiked_' name ext];
        for v=1:size(Y,4)
            V(v).descrip = 'spmup despiked';
            V(v).fname=fname;
            despiked = spm_write_vol(V(v),squeeze(YY(:,:,:,v)));
        end
    end
end

%% write the report for QA
disp('saving spmup_despike_log')
if ischar(P); spmup_despike_log.P = P; end
spmup_despike_log.flags = flags;
if isempty(flags.window)
    spmup_despike_log.window = R;
else
    spmup_despike_log.window = window;
end
Despiked_QA = sqrt(mean((YY-Y).^2,4)); % how much different
spmup_despike_log.RMS = Despiked_QA;
[pathstr,name,ext]= fileparts(V(1).fname);
save([pathstr filesep 'spmup_despike_log'],'spmup_despike_log')

% figure('Name','QA'); colormap('hot');
% for z=1:size(Despiked_QA,3)
%     imagesc(flipud(Despiked_QA(:,:,z)')); axis square;
%     title(['Slice ' num2str(z)]); pause
% end

disp('despiking done')
disp('--------------')