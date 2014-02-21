function spmup_despike(varargin)
% 
% SPM U+ routine to 'despike' fMRI time-series in a similar way as AFNI does
% Note is requires the statistics toolbox (nansum, icdf are called)
%
% FORMAT spmup_despike
%        spmup_despike(P)
%        spmup_despike(P,M)
%        spmup_despike(P,M,flags)
%        spmup_despike(P,[],flags)
%
% INPUT if none the user is prompted
%       P the names of the fMRI images (time-series)
%       M the name of the mask
%       flags defines options to be used
%             'flags.auto_mask','off' or 'on' if M is not provided, auto_mask is
%             'on' but if set to 'off' the user is prompted to select a mask
%
% OUTPUT spmup_despike_log is saved onto disk where the data are
%        spmup_despike_log can be reviewed using spmup_review_despike_log
%        spmup_despike_log is structure with the fields:
%        - P the list of files (data) analyzed
%        - outlying_voxels the proportion of outlying voxels per volume
%        - outlying_volumes a binary classification of volumes
%        - despiked_voxels the proportion of despiked voxels per volume
%        - class a 4d binary matrix indicating despiked voxels
%
%        the new dataset with the spikes removed is also written in the
%        data folder with the prefix 'despiked_'
%
% --------------------------------------------------------------------------
% First, one look for outlier voxels based on the Median Absolute Deviation
% Here the MAD is median absolute value of time series minus trend.
% The trend is optained using son_detrend (2nd order polynomial)
% Outliers are then those voxels with values outside 
% k = alphav * sqrt(pi/2) * MAD  --- this is the similar as 3dToutcount
%
% Second, if some volume outliers are found, despike
% Note that the despiking is not based on the results of the above, some 
% voxels from volume not seen as outliers can still be despiked.
% The data are 1st smoothed either with a median filter (N/30 points window)
% or using a smooth curve (as in AFNI - see flag options) and then we look
% for outliers in the residuals and the data are interpolated as in AFNI.
% Finaly, we write the data with the prefix despiked_  
% --------------------------------------------------------------------------
%
% Cyril Pernet 17 Feb 2014
% --------------------------------------------------------------------------
% Copyright (c) SPM U+ toolbox

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
        if strcmp(varargin(3).auto_mask,'on') || strcmp(varargin(3).auto_mask,'off')
            flags.auto_mask = varargin(3).auto_mask;
        else
            error('flags.auto_mask must be ''on'' or ''off''');
        end
    end
    
    if isfield(varargin{3},'method')
        if strcmp(varargin(3).method,'median') || strcmp(varargin(3).method,'curve')
            flags.method = varargin(3).method;
        else
            error('flags.method must be ''median'' or ''curve''');
        end
        
        if strcmp(varargin(3).method,'curve')
            error('oops curve fitting is not implemented yet')
         end
    end
end


%% get data and mask 
% memory mapped data
if get_data == 1;
    [P,sts] = spm_select(Inf,'image','select the time seires',[],pwd,'.*',1);
else
    P = varargin{1};
end
V = spm_vol(P);

% memory mapped mask
if get_mask == 1;
    [M,sts] = spm_select(1,'image','select the mask',[],pwd,'.*',1);
    if sts ~=0
       Mask = spm_read_vols(spm_vol(M));
    else
        error('spm u+ stopped - mask selection interupted')
    end
elseif get_mask == 0 && strcmp(flags.auto_mask,'off')
    try
        Mask = spm_read_vols(spm_vol(varargin{2}));
    catch varin_error
        error(varin_error.name)
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

% data
Y = spm_read_vols(V);
N = size(Y,4);

%% this part is similar as 3dToutcount 
% http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dToutcount.html

% detrend, get the MAD and classify
class = NaN(size(Y));
alphav = icdf('Normal',1-(0.001/N),0,1);

disp('checking outliers')
for x=1:size(Y,1)
    for y=1:size(Y,2)
        index = find(Mask(x,y,:));
        if ~isempty(index) 
        clean_data = spm_detrend(squeeze(Y(x,y,:,:)),2); % detrend
        M = repmat(nanmedian(clean_data,2),1,N); % medians of the time-series
        MAD = median(abs(clean_data-M),2); % Median absolute deviation of the time series
        k = (alphav * sqrt(pi/2)) .* MAD; % how far is far away
        up = repmat(mean(clean_data,2)+k,1,N);
        down = repmat(mean(clean_data,2)-k,1,N);
        class(x,y,:,:) = (clean_data > up) + (clean_data < down); % threshold
        end
    end
end
 
% compute the proportion of outliers per volume
Nb_voxels = nansum(Mask(:));
for im=1:N
    tmp = squeeze(class(:,:,:,im));
    outlying_voxels(im) = (nansum(tmp(:))./Nb_voxels)*100;
end
M = repmat(median(outlying_voxels),1,N);
MAD = median(abs(outlying_voxels-M));
MADN = repmat((MAD./.6745),1,N); % this is almost like 3.5*MAD but better
outlying_volumes = (abs(outlying_voxels-M) ./ MADN) > sqrt(chi2inv(0.975,1));

%% now do the despiking
% this part is similar to 3dDespike execpt we use a simple median smoother
% http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html

YY =Y;
class = NaN(size(Y));
despiked_voxels = [];
if sum(outlying_volumes) ~=0
    index = 1; disp('despiking data ... '); 
    tot = size(Y,1)*size(Y,2)*size(Y,3);
    f = waitbar(0,'Percentage done','name','Despiking');
    for x=1:size(Y,1)
        for y=1:size(Y,2)
            for z=1:size(Y,3) % smooth the data voxel wise
                % fprintf('percentafge done: %g \n',index/tot*100) 
                waitbar((index)/tot); index = index+1;
                if Mask(x,y,z) ~= 0 && ~isnan(Mask(x,y,z))
                    data = squeeze(Y(x,y,z,:))';
                    newdata = nan(size(data));
                    
                    if mod(ceil(N/30),2) % % as in AFNI /30
                        window = ceil(N/30);
                    else
                        window = floor(N/30);
                    end
                    
                    if window < 3 % need at least 3 points
                        window = 3;
                    end
                    
                    % median smoothing
                    % ------------------
                    if strcmp(flags.method,'median')
                        % beginning
                        for p=1:floor(window/2)
                            newdata(p) = median([repmat(data(1),1,ceil(window/2)-p) data(1:p) data(p+1:p+floor(window/2))]);
                        end
                        % middle
                        for p=(floor(window/2)+1):(N-floor(window/2))
                            newdata(p) = median(data((p-floor(window/2)):(p+floor(window/2))));
                        end
                        % end
                        last = 1;
                        for p=(N-floor(window/2)+1):(N-1) % don't do last data point
                            newdata(p) = median([data(p-ceil(window/2):p-1) repmat(data(p),1,floor(window/2)-last)]);
                            last = last+1;
                        end
                        newdata(N) = data(N);
                        
                    % N/30 points curve fitting (1st order)
                    % ------------------------------------
                    elseif strcmp(flags.method,'curve')
                        
                        % f(t) = a+b*t+c*t^2+cumsum(d*sin((2*pi*t)/T)+e*cos((2*pi*t)/T))
                        % T duration of the time series
                        % a,b,c,d,e are minimized for v(t)-f(t)
                        % TO DO - use fminsearch on that function
                        
                    end
                    
                    
                    % MAD of the residuals
                    res = data-newdata;
                    MAD = median(abs(res - repmat(median(res),1,N)));
                    SIGMA = sqrt(pi/2)*MAD;
                    s = res/SIGMA;
                    
                    %  * Values with s > c1 are replaced with a value that yields
                    %     a modified s' = c1+(c2-c1)*tanh((s-c1)/(c2-c1)).
                    %  * c1 is the threshold value of s for a 'spike' [default c1=2.5].
                    %  * c2 is the upper range of the allowed deviation from the curve:
                    %     s=[c1..infinity) is mapped to s'=[c1..c2)   [default c2=4].
                    out = find(s > 2.5);
                    class(x,y,z,out) = 1;
                    c1 = 2.5; c2=4; s2 = s;
                    for p=1:length(out)
                        s2(out(p)) = c1+(c2-c1)*tanh((s(out(p))-c1)/(c2-c1));
                    end
                    
                    % reverse s2 to the real data
                    YY(x,y,z,:) = (s2*SIGMA)+newdata;
                end
            end
        end
    end
    close(f)
    
    for im=1:N
        tmp = squeeze(class(:,:,:,im));
        despiked_voxels(im) = (nansum(tmp(:))./Nb_voxels)*100;
    end
    
    figure('Name','Despiking')
    subplot(2,2,[1 2]);
    plot(outlying_voxels,'LineWidth',3); grid on;
    xlabel('volumes','Fontsize',12); axis tight
    ylabel('percentage of outlying voxels','Fontsize',12);
    title('Volume Outlier detection','Fontsize',14); hold on
    plot(outlying_volumes.*(max(outlying_voxels)-mean(outlying_voxels))+mean(outlying_voxels),'or','LineWidth',2);
    subplot(2,2,[3 4]);
    plot(despiked_voxels,'LineWidth',3); grid on;
    xlabel('volumes','Fontsize',12); axis tight
    ylabel('percentage of despiked voxels ','Fontsize',12);
    title('Voxel despiked','Fontsize',14); drawnow
    saveas(gcf, 'despiking - volume outlier detection.eps','psc2'); close(gcf)
else
    disp('no outlying fund, data are saved like the originals .. ')
    figure('Name','Proportion of outlying voxels')
    plot(outlying_voxels,'LineWidth',3); grid on;
    xlabel('volumes','Fontsize',12); axis tight
    ylabel('percentage of outlying voxels','Fontsize',12);
    title('Volume Outlier detection - No despiking performed','Fontsize',14); 
    drawnow; saveas(gcf, 'despiking - volume outlier detection.eps','psc2'); close(gcf)
end

%% write the data
for v=1:size(Y,4)
    disp('writing data')
    V(v).descrip = 'spmu+ despiked';
    [pathstr,name,ext]= fileparts(V(v).fname);
    V(v).fname = [pathstr filesep 'despiked_' name ext];
    spm_write_vol(V(v),squeeze(YY(:,:,:,v)));
end

%% write the report
disp('saving spmup_despike_log')
spmup_despike_log.P = P;
spmup_despike_log.outlying_voxels = outlying_voxels;
spmup_despike_log.outlying_volumes = outlying_volumes;
spmup_despike_log.despiked_voxels = despiked_voxels;
spmup_despike_log.class = class;
save spmup_despike_log spmup_despike_log
disp('despiking done')



