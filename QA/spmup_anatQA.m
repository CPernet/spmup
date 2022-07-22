function [anatQA,int_data] = spmup_anatQA(varargin)

% *This implements some of the the anatomical QA measures from the*
% *Preprocessed Connectome Project Quality Assurance Protocol (QAP)*
% <http://preprocessed-connectomes-project.org/quality-assessment-protocol/>
% 
% FORMAT: anatQA = spmup_anatQA(anat,c1,c2,option)
% 
% INPUT: anat is the anatomical image
%        c1 is the gray matter image
%        c2 is the white matter image
%        option is 'no background' set to 'true' or 'false'
%                  'figure': off/on/save (default) 
%
% OUTPUT anatQA is a structure with the following fields:
%        - SNR : the signal-to-Noise Ratio, ie the mean intensity within gray and white matter divided
%                by the standard deviation of the values outside the brain*. Higher values are better.
%        - CNR : the Contrast to Noise Ratio, i.e. the mean of the white matter intensity values minus the mean 
%                of the gray matter intensity values divided by the standard deviation of the values outside the
%                brain*. Higher values are better.
%        - FBER: Foreground to Background Energy Ratio, i.e. the variance of voxels in grey and white matter 
%                divided by the variance of voxels outside the brain*. Higher values are better.
%
%        - EFC : Entropy Focus Criterion, i.e. the entropy of voxel intensities proportional to the maximum 
%                possible entropy for a similarly sized image. Indicates ghosting and head motion-induced blurring. 
%                Lower values are better. See <http://ieeexplore.ieee.org/document/650886/>
%
% * when the background has 0 variance (e.g. a sequence with noise suppression like FLAIR) then the standard 
%   deviation of the white matter is used as reference instead of the background
%
% Note: GM and WM are thresholded by making them mutually exclusive
% The background is found using data outside a large brain mask and
% trimming extreme values
% 
% Cyril Pernet 20 January 2017
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2017

% 2nd output int_data is used for unit testing, this retuns in a structure
% intermediate values   

%% input data
if nargin == 0
    [Anat,sts] = spm_select(1,'image' ,'Select Anatomical Image',{},pwd,'.*','1');
    if sts == 1
        AnatV  = spm_vol(Anat);
    else
        return
    end
    
    [Gray,sts] = spm_select(1,'image' ,'Select Gray matter image',{},pwd,'c1.*','1');
    if sts == 1
        GrayV  = spm_vol(Gray);
    else
        return
    end
    
    [White,sts] = spm_select(1,'image' ,'Select White matter image',{},pwd,'c2.*','1');
    if sts == 1
        WhiteV  = spm_vol(White);
    else
        return
    end
else
    AnatV  = spm_vol(varargin{1});
    GrayV  = spm_vol(varargin{2});
    WhiteV = spm_vol(varargin{3});
end

if any(AnatV.dim ~= GrayV.dim)
    error('Ther anatomical and gray matter images do not have the same size')
end

if any(AnatV.dim ~= WhiteV.dim)
    error('Ther anatomical and white matter images do not have the same size')
end

option.background = 1;
option.fig        = 'save'; % make invisible figures and save as/append spmup_QC.ps
if nargin > 3
    for n = 4:2:nargin
        if strcmpi(varargin{n},'no background') && strcmpi(varargin{n+1},'true')
            option.background = 0;
        elseif contains(varargin{n},'fig')
            option.fig  = varargin{n+1}; 
        end
    end
end
    
    
%% compute

disp('---------------------------')
fprintf(' running spmup_anatQA on %s\n',AnatV.fname)
disp('---------------------------')

% define a large brain mask making sure the background has nothing but background in it
brain_mask = (smooth3(spm_read_vols(GrayV),'box',25)+smooth3(spm_read_vols(WhiteV),'box',25)) > 0.1;
[x,y,z]    = ind2sub(AnatV.dim,find(brain_mask==0));
data       = sort(spm_get_data(AnatV,[x y z]'));
if isempty(data)
    figure; imagesc(squeeze(brain_mask(:,:,round(AnatV.dim(3)/2)))); title('brain mask')
    error('no background data were obtained from the brain mask \n error in image %s\n',AnatV.fname);
end

% compute the noise level using and 50% trimmed data (remove extreme values)
data(isnan(data)) = [];
up                = median(data)+iqr(data)/2; 
down              = median(data)-iqr(data)/2;
index             = logical((data<up).*(data>down));
std_nonbrain      = std(data(index)); 
if std_nonbrain == 0
    option.background = 0;
end

%% make a figure showing the mask and the data range
if ~strcmpi(option.fig,'off')
    figure('Name','Brain Mask')
    if strcmpi(option.fig,'on')
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    else
        set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1],'visible','off')
    end
    
    subplot(2,2,1);
    tmp                = spm_read_vols(AnatV); % read anat
    tmp                = tmp ./max(tmp(:))*100; % nornalize to 100
    tmp(brain_mask==0) = 200; % set non brain to 200
    imagesc(fliplr(squeeze(tmp(:,:,round(AnatV.dim(3)/2))))');
    title('brain mask - sagital');
    subplot(2,2,2);
    imagesc(flipud(squeeze(tmp(:,round(AnatV.dim(2)/2),:))'));
    title('brain mask - axial');
    subplot(2,2,4);
    imagesc(flipud(squeeze(tmp(round(AnatV.dim(1)/2),:,:))'));
    title('brain mask - coronal');
    clear tmp
    
    subplot(2,2,3);
    k = round(1+log2(length(data)));
    histogram(data,k,'FaceColor',[1 1 0],'EdgeColor',[1 1 0],'FaceAlpha',0.8); hold on
    plot(repmat(down,max(hist(data,k)),1),[1:max(hist(data,k))],'k','LineWidth',3);
    plot(repmat(up,  max(hist(data,k)),1),[1:max(hist(data,k))],'k','LineWidth',3);
    clear x y z; [x,y,z] = ind2sub(AnatV.dim,find(brain_mask));
    histogram(spm_get_data(AnatV,[x y z]'),k,'FaceColor',[0 0 1],'EdgeColor',[0 0 1],'FaceAlpha',0.4); hold on
    title('background voxels & limits vs brain mask'); axis tight; box on; grid on;
    
    [filepath,filename] = fileparts(AnatV.fname);
    if isempty(filepath)
        filepath = pwd;
    end
    
    if exist(fullfile(filepath,'spmup_QC.ps'),'file')
        print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(filepath,'spmup_QC.ps'));
    else
        print (gcf,'-dpsc2', '-bestfit', fullfile(filepath,'spmup_QC.ps'));
    end
end

%% now do all computations
% make gray/white mutually exclusive
I_GM = spm_read_vols(GrayV)  > spm_read_vols(WhiteV);
I_WM = spm_read_vols(WhiteV) > spm_read_vols(GrayV); 

% compute taking voxels in I_GM and I_WM
clear x y z 
[x,y,z] = ind2sub(AnatV.dim,find(I_GM));
dataGM  = spm_get_data(AnatV,[x y z]');
meanGM  = mean(dataGM);

clear x y z 
[x,y,z] = ind2sub(AnatV.dim,find(I_WM));
dataWM  = spm_get_data(AnatV,[x y z]');
meanWM  = mean(dataWM);

if option.background == 0
    anatQA.SNR  = ((meanGM+meanWM)/2)/std(dataWM);
    anatQA.CNR  = (meanWM-meanGM)/std(dataWM);
    anatQA.FBER = var([dataGM dataWM])/std(dataWM)^2;
else
    anatQA.SNR  = ((meanGM+meanWM)/2)/std_nonbrain;
    anatQA.CNR  = (meanWM-meanGM)/std_nonbrain;
    anatQA.FBER = var([dataGM dataWM])/std_nonbrain^2;
end

data = spm_read_vols(AnatV);
data(isnan(data)) = [];
Bmax = sqrt(sum(data(:).^2));
try
    % should work most of the time, possibly throwing warnings (cf Joost Kuijer)
    anatQA.EFC = real(nansum((data(:)./Bmax).*log((data(:)./Bmax))));
catch
    anatQA.EFC = real(nansum((data(:)./Bmax).*abs(log((data(:)./Bmax)))));
end

%% outpouts for unit test
if nargout == 2
    int_data.meanGM = meanGM;
    int_data.meanWM = meanWM;
    int_data.std_nonbrain = std_nonbrain;
    int_data.varGMWM = var([dataGM dataWM]);
    int_data.Bmax = Bmax;
end

if nargout == 0
    writetable(struct2table(anatQA), fullfile(filepath,[filename '_anatQA.txt']));
end

