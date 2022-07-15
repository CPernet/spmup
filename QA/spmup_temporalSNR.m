function tSNR = spmup_temporalSNR(time_series,masks,varargin)

% Computes the temporal SNR of the time_series input in the different
% compartments provided by the masks images = mean signal / std over time
% The routine recapitulates tSNR as described in Thomas Liu  (2016),
% Caesar Caballero-Gaudes and Richard C. Reynolds (2016). &
% Lawrence Wald and Jonathan R Polimeni (2016).
%
% FORMAT tSNR = spmup_temporalSNR(time_series,masks,options)
%
% USAGE: tSNR = spmup_temporalSNR('4Dimg.nii',{'GM.nii','WM.nii','CSF.nii','Mask.nii'},...
%                   'figure','save','SNR0','on','linearity','off')
%
% INPUTS - time_series: a cell array of file names (see spm_select)
%       - masks: a cell array of (non binary) tissue class images (GM, WM,
%       CSF in that order) in the same space as the time series (i.e.
%       typically from the anatomical coregistered to the mean EPI) +
%       optional a brain mask (if not provided, one is computed as the sum 
%       of smoothed tissues > 0.5).
%       - Options
%       'figure': off/on/save (default) 
%       'SNR0': off (default)/on
%       'linearity': off (default)/on
%
% OUTPUTS - tSNR is a structure with the following fields:
%            .GM: mean GM signal / std over time (estimate BOLD from GM>(WM+CSF))
%            .WM:  mean WM signal / std over time (estimate non-BOLD from WM>(GM+CSF))
%            .CSF: mean CSF signal / std over time (estimate non-BOLD from CSF>(GM+WM))
%            .Background:  mean signal outside mask (GM+WM+CSF) / std over time
%            .average (tSNR): mean signal / sqrt(std(GM)^2+std(WM+CSF)^2+std(Background)^2)
%             OPTIONAL
%            .image (SNR0):  mean signal inside mask / std outside mask over time
%            .physio2termal_ratio: sqrt((tSNR(whole image)/SNR0(brain only))^2-1)
%            .physio2termal_corr: correlation between images
%            .signal_mean: sqrt(std(GM)^2+std(WM+CSF)^2) / sqrt((SNR0^2/tSNR- 1)/SNR0^2)
%                          Since tSNR = SNR0^2 / (1+L^2*SNR0^2), 
%                          we have L^2 = (SNR0^2 /tSNR - 1) / SNR0^2
%                          and sqrt(std(GM)^2+std(WM+CSF)^2) = L*Smean
%            .roi: tSNR for increased ROI (from in mask by increasing slices) ~linear function of srqrt(nb voxels)
%            Since tSNR = SNR0^2 / (1+L^2*SNR0^2), we have L^2 = (SNR0^2 /tSNR - 1) / SNR0^2
%            and sqrt(std(GM)^2+std(WM+CSF)^2) = L*Smean
%        - tSNR_time_series.nii image is also saved on the drive, showing tSNR in each voxel for GM, WM and CSF as computed above
%
%        tSNR_time_series.nii image is also saved on the drive, showing tSNR in each voxel 
%        for GM, WM and CSF as computed above
%
% References: (1) Thomas Liu  (2016) 
%                 Noise contributions to the fMRI signal: An overview
%                 NeuroImage, 143, 141-151.
%             (2) Cesar Caballero-Gaudes and Richard C. Reynolds (2016). 
%                 Methods For Cleaning The BOLD fMRI Signal. 
%                 NeuroImage, 154,128-149
%             (3) Lawrence Wald and Jonathan R Polimeni (2016). 
%                 Impacting the effect of fMRI noise through hardware and acquisition choices – 
%                 Implications for controlling false positive rates. 
%                 NeuroImage, 154,15-22
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

%% defaults
fig  = 'save'; % make invisible figures and save as/append spmup_QC.ps 
snr0 = 'on';   % computes added estimates signal vs thermal noise
roi  = 'off';  % do not compute the linearily of tSNR vs roi size (takes time)

%% check inputs
if nargin == 0
    [time_series,sts] = spm_select(Inf,'image' ,'Select your fMRI time series',{},pwd,'.*',Inf);
    if sts == 0
        return
    end
    [masks,sts]       = spm_select(3,'image' ,'Select tissue masks (GM,WM,CSF)',{},pwd,'.*',1);
    if sts == 0
        return
    end
end

if size(time_series,1) == 1 && strcmp(time_series(length(time_series)-1:end),',1')
    time_series = time_series(1:length(time_series)-2); % in case picked 4D put left ,1
end

if iscell(time_series)
    for v=size(time_series,1):-1:1
        V(v) =spm_vol(time_series{v});
    end
else
    V = spm_vol(time_series);
end

if iscell(V); V = cell2mat(V); end
if size(V,1) < 10; error('there is less than 10 images in your time series ??'); end

if iscell(masks) && size(masks,2) > size(masks,1)
    masks = masks'; % appears like char array P but as cell
end
    
for m=size(masks,1):-1:1
    if iscell(masks) ==0
        VM(m) = spm_vol(masks(m,:));
    else
        VM(m) = spm_vol(masks{m});
    end
end

if iscell(VM); VM = cell2mat(VM); end

if length(VM) < 3
    error(['at least 3 masks files expected, ' num2str(size(VM,1)) ' detected - check input file'])
end

if any(VM(1).dim~= VM(2).dim) || any(VM(3).dim~= VM(2).dim)
    error('masks do not have the same dimentions')
end

if any(V(1).dim~= VM(1).dim)
    error('dimention mismatch between the time series and masks')
end

if nargin > 2
   for o = 1:length(varargin)
      if strcmpi(varargin{o},'figure')
          fig = varargin{o+1};
      elseif strcmpi(varargin{o},'SNR0')
          snr0 = varargin{o+1};
      elseif strcmpi(varargin{o},'linearity')
          roi = varargin{o+1};
      end
   end
end

[~,filename] = fileparts(V(1).fname);
disp   ('-------------------------------')
fprintf(' running spmup_temporalSNR on %s\n',filename)
disp   ('-------------------------------')

%% Compute relative masks
GM         = spm_read_vols(VM(1));
WM         = spm_read_vols(VM(2));
CSF        = spm_read_vols(VM(3));
if length(VM) == 4
    brain_mask = spm_read_vols(VM(4));
else
    brain_mask = (smooth3(GM,'box',15)+smooth3(WM,'box',15)+smooth3(CSF,'box',15))>0;
end
GM         = GM.*(GM>0.5);
WM         = WM.*(WM>0.5);
CSF        = CSF.*(CSF>0.5);     % baseline prob 50%
GM         = GM.*(GM>(WM+CSF));  % figure; rst_hist(GM(:))
WM         = WM.*(WM>(GM+CSF));  % figure; rst_hist(WM(:))
CSF        = CSF.*(CSF>(WM+GM)); % figure; rst_hist(CSF(:))

%% in masks tSNR
clear x y z
[x,y,z]  = ind2sub(size(GM),find(GM));
data     = spm_get_data(V,[x y z]');
stdGM    = nanmean(nanstd(data,1));
tSNR.GM  = nanmean(nanmean(data,1)) /stdGM;

clear x y z
[x,y,z]  = ind2sub(size(WM),find(WM));
data     = spm_get_data(V,[x y z]');
stdWM    = nanmean(nanstd(data,1));
tSNR.WM  = nanmean(nanmean(data,1)) /stdWM;

clear x y z
[x,y,z]  = ind2sub(size(CSF),find(CSF));
data     = spm_get_data(V,[x y z]');
stdCSF   = nanmean(nanstd(data,1));
tSNR.CSF = nanmean(nanmean(data,1)) /stdCSF;

%% Background
clear x y z
[x,y,z]             = ind2sub(size(brain_mask),find(brain_mask ~= 1));
data                = spm_get_data(V,[x y z]');
if sum(data(:)) ~= 0
    try
        stdBackground       = nanmean(nanstd(data,1));
    catch
        warning('serious issue can''t get std in background! trying voxel-wise')
        out = zeros(1,size(data,2));
        for v=1:size(data,2)
            try
                out(v) = nanstd(squeeze(data(:,v)),1);
            catch
                fprintf('serious issue voxel %g %g %g \n',x(v),y(v),z(v))
            end
        end
        stdBackground       = nanmean(out);
    end
else
    % data have been masked - backgound = 0, looking for a mean image
    % might not work either
    [fpath,fname,ext] = fileparts(V(1).fname);
    if exist(fullfile(fpath,['mean' fname ext]),'file')
        data          = spm_get_data(spm_vol(fullfile(fpath,['mean' fname ext])),[x y z]');
        stdBackground = nanmean(nanstd(data,1));
    else
        stdBackground = 0;
    end
end

% tSNR.Background_raw = data;
tSNR.Background     = nanmean(nanmean(data,1)) /stdBackground; % figure; rst_hist(data(:))

% Computes the density estimate of data using a Random Average Shifted
% Histogram algorithm is coded based on Bourel et al. Computational
% Statistics and Data Analysis 79 (2014)
[filepath, filename, ext] = fileparts(time_series(1,:));
if contains(ext,'nii') % make sure we remove any frame trail
    ext = '.nii';
elseif contains(ext,'img')
    ext = '.img';
end

if sum(isnan(data(:))) ~= numel(data) && stdBackground ~= 0 && ...
        ~strcmpi(fig,'off')
         
    % tSNR per voxel from background
    data = (nanmean(data,1)/stdBackground)';
    
    % see where is it and spatial distribition
    SNRimage = zeros(V(1).dim); index = 1;
    SNRimage(find(brain_mask~=1)) = data;
    mymin = nanmedian(data)-3*iqr(data);
    if mymin<0 || isnan(mymin)
        mymin = 0;
    end

    figure('Name','Background SNR');
    if strcmpi(fig,'on')
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    else
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1],'visible','off')
    end
    figindex = [1 2 3 4 10 11 12 13 19 20 21 22 28 29 30 31];
    for z=1:floor(V(1).dim(3)./16)+1:V(1).dim(3)-1
        subplot(4,9,figindex(index));
        imagesc((squeeze(SNRimage(:,:,z))'));
        index = index+1; colormap(cubehelix(32,[1.15,0.1,4,1], [0,1], [0,0.85]))
        try caxis([mymin, nanmedian(data)+3*iqr(data)]); end
        set(gca,'XtickLabel',[],'YtickLabel',[])
        xlabel(['slice ' num2str(z)]);
    end
    clear SNRimage
    
    % do the histogram (kernel density)
    n     = length(data);
    m     = 100; % number of hist to compute
    h     = 2.15*sqrt(var(data))*n^(-1/5);
    delta = h/m;
    % 1 make a mesh with size delta
    t0      = min(data)-(min(diff(data))/2);
    tf      = max(data)+(min(diff(data))/2);
    nbin    = ceil((tf-t0)/delta);
    binedge = t0:delta:(t0+delta*nbin);
    out     = find(binedge>tf);
    if isempty(out) || out == 1
        binedge(out) = tf;
    else
        binedge(out(1))     = tf;
        binedge(out(2:end)) = [];
    end
    % 2 Get the weight vector.
    kern = @(x) (15/16)*(1-x.^2).^2; % inline('(15/16)*(1-x.^2).^2');
    ind  = (1-m):(m-1);
    den  = sum(kern(ind/m));% Get the denominator.
    wm   = m*(kern(ind/m))/den;% Create the weight vector.
    % 3 compute bin with shifted edges
    RH  = zeros(1,nbin);
    RSH = zeros(m,nbin);
    for e=1:m
        v        = binedge + (delta*randn(1,1)); % e is taken from N(0,h);
        v(v<t0)  = t0; % lower bound
        v(v>tf)  = tf; % upper bound
        nu       = histc(data,v);
        nu       = [zeros(1,m-1) nu' zeros(1,m-1)];
        for k = 1:nbin
            ind  = k:(2*m+k-2);
            RH(k)= sum(wm.*nu(ind));
        end
        RSH(e,:) = RH/(n*h);
    end
    K  = mean(RSH,1);
    bc = t0+((1:nbin)-0.5)*delta;
    
    subplot(4,8,[5 6 7 8 13 14 15 16 21 22 23 24 29 30 31 32]);
    bar(bc,K,1,'FaceColor',[0.5 0.5 1]);
    title('RAS Histogram - background noise');
    grid on; box on; ylabel('tSNR'); drawnow
    if exist(fullfile(filepath,'spmup_QC.ps'),'file')
        print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(filepath,'spmup_QC.ps'));
    else
        print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(filepath,'spmup_QC.ps'));
    end
    close(gcf)
end

%% average (tSNR)
[x,y,z]  = ind2sub(size(WM),find(WM+CSF));
data     = spm_get_data(V,[x y z]');
stdWMCSF = nanmean(nanstd(data,1)); % presumably non BOLD
[x,y,z]  = ind2sub(size(brain_mask),find(GM+WM+CSF+(brain_mask ~= 1)));
try
    data         = spm_get_data(V,[x y z]'); % the whole image or so
    tSNR.average = nanmean(nanmean(data,1)) / sqrt(stdGM^2+stdWMCSF^2+stdBackground^2);
    data         = (nanmean(data,1)/sqrt(stdGM^2+stdWMCSF^2))';
catch
    warning('serious issue can''t get all data! trying voxel-wise')
    data = zeros(V(1).dim);
    for v=length(x):-1:1
        try
            tmp                  = spm_get_data(V,[x(v) y(v) z(v)]');
            average(v)           = nanmean(tmp);
            data(x(v),y(v),z(v)) = nanmean(tmp)/sqrt(stdGM^2+stdWMCSF^2);
        catch
            fprintf('serious issue voxel %g %g %g \n',x(v),y(v),z(v));
            average(v) = NaN;
            data(x(v),y(v),z(v)) = 0;
        end
    end
    tSNR.average = nanmean(average) / sqrt(stdGM^2+stdWMCSF^2+stdBackground^2);
end
SNRimage     = zeros(V(1).dim);
SNRimage(find(GM+WM+CSF+(brain_mask ~= 1))) = data;

W                    = V(1);
W.fname              = [filepath filesep filename(1:end-5) '_tSNR' ext];
W.private.dat.fname  = [filepath filesep filename(1:end-5) '_tSNR' ext];
W.descrip            = 'tSNR image - see spmup_temporalSNR';
W.private.dat.dim    = V(1).private.dat.dim(1:3);
W.n                  = [1 1];
spm_write_vol(W,SNRimage);

%% --------  optional -------------------

if strcmpi(snr0,'on') && stdBackground ~= 0
    %% SNR0 (brain only)
    [x,y,z]    = ind2sub(size(brain_mask),find(GM+WM+CSF));
    data       = spm_get_data(V,[x y z]');
    tSNR.image = nanmean(nanmean(data,1)) / stdBackground;
    data       = (nanmean(data,1)/stdBackground)';
    SNROimage  = zeros(V(1).dim);
    SNROimage(find(GM+WM+CSF)) = data;
    
    if ~strcmpi(fig,'off')
      
        if strcmpi(fig,'on')
            figure('Name','SNR0');
        else
            figure('Name','SNR0','Visible','Off');
        end
        figindex = [1 2 3 4 10 11 12 13 19 20 21 22 28 29 30 31];
        index = 1;
        for z=1:floor(V(1).dim(3)./16)+1:V(1).dim(3)-1
            subplot(4,9,figindex(index));
            imagesc(flipud(squeeze(SNRimage(:,:,z)')));
            caxis([min(SNRimage(:)), max(SNRimage(:))]); set(gca,'XtickLabel',[],'YtickLabel',[]);
            xlabel(['slice ' num2str(z)]); if index ==2; title('tSNR image'); end
            subplot(4,9,figindex(index)+5); R= abs(sqrt((squeeze(SNROimage(:,:,z)'./SNRimage(:,:,z)').^2)-1));
            imagesc(flipud(R)); if index ==2; title('SNR0/tSNR image'); end
            try caxis([min(R(:)), max(R(:))]); end
            set(gca,'XtickLabel',[],'YtickLabel',[]); xlabel(['slice ' num2str(z)]);
            index = index+1; colormap(cubehelix(32,[3,1.9,1.5,1], [0,1], [0.2,0.8]))
        end
        drawnow
        if exist(fullfile(filepath,'spmup_QC.ps'),'file')
            print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(filepath,'spmup_QC.ps'));
        else
            print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(filepath,'spmup_QC.ps'));
        end
        close(gcf)
    end
    
    tSNR.SNR02tSNR_corr = corr(SNRimage(:),SNROimage(:));
    clear SNRimage SNROimage
    
    %% ratio
    tSNR.SNR02tSNR_ratio = abs(sqrt(((tSNR.average/tSNR.image)^2)-1));
    
    %% signal
    L2 = (tSNR.image^2 /tSNR.average - 1) / tSNR.image^2;
    tSNR.signal_mean = sqrt(stdGM^2+stdWMCSF^2) / sqrt(L2);
end

if strcmpi(roi,'on')
    %% per ROI (absolute masking)
    disp('tSNR - checking linearity ..')
    GM  = spm_read_vols(VM(1));
    WM  = spm_read_vols(VM(2));
    CSF = spm_read_vols(VM(3));
    
    index = 1;
    for p=0.95:-0.05:0.1
        
        [x,y,z] = ind2sub(size(GM),find(GM>=p));
        data    = spm_get_data(V,[x y z]');
        stdGM   = nanmean(nanstd(data,1));
        
        [x,y,z]  = ind2sub(size(WM),find(WM>p+CSF>p));
        data     = spm_get_data(V,[x y z]');
        stdWMCSF = nanmean(nanstd(data,1));
        
        ROI = (GM>p+WM>p+CSF>p);
        % figure; for z=1:size(ROI,3); imagesc(squeeze(ROI(:,:,z))); pause; end
        [x,y,z]               = ind2sub(size(brain_mask),find(ROI+(brain_mask~=1)));
        data                  = spm_get_data(V,[x y z]'); % the whole image or so
        tSNR.roi.value(index) = nanmean(nanmean(data,1)) / sqrt(stdGM^2+stdWMCSF^2+stdBackground^2);
        tSNR.roi.size(index)  = nansum(ROI(:));
        index = index + 1;
    end
    % fit a line to this
    B = pinv([sqrt(tSNR.roi.size)' ones(18,1)])*tSNR.roi.value';
    model = [sqrt(tSNR.roi.size)' ones(18,1)]*B;
    tSNR.roi.slope = B(1);
    
    if ~strcmpi(fig,'off')
      
        if strcmpi(fig,'on')
            figure('Name','SNR per size');
        else
            figure('Name','SNR per size','Visible','Off');
        end
        plot(sqrt(tSNR.roi.size),[sqrt(tSNR.roi.size)' ones(18,1)]*B,'LineWidth',3);
        hold on; plot(sqrt(tSNR.roi.size),tSNR.roi.value,'ro','LineWidth',2);
        axis tight; box on; grid minor; ylabel('temporal SNR','FontSize',12)
        xlabel('sqrt of the number of in brain voxels used','FontSize',12)
        mytitle = sprintf('tSNR=%g*sqrt(nb of voxels)+%g \n RMSE=%g',B(1),B(2),sqrt(mean(model - tSNR.roi.value')));
        title(mytitle,'FontSize',12); drawnow
        if strcmpi(fig,'save')
            if exist(fullfile(filepath,'spmup_QC.ps'),'file')
                print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(filepath,'spmup_QC.ps'));
            else
                print (gcf,'-dpsc2', '-bestfit', fullfile(filepath,'spmup_QC.ps'));
            end
            close(gcf)
        end
    end
end

