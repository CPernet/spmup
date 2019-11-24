function [FD,RMS,FD_outliers,RMS_outliers] = spmup_FD(varargin)

% simple routine that computes framewise displacement as a sum (FD) and RMS 
% Power et al. (2012) doi:10.1016/j.neuroimage.2011.10.018  
% Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048
%
% FORMAT: [FD,RMS,motion] = spmup_FD
%         [FD,RMS,motion] = spmup_FD(realignment_file)
%         [FD,RMS,motion] = spmup_FD(realignment_file,'Radius',value,'figure','on/save/off')
%
% INPUT: realigment is the realigment parameters file
%        radius is the expected head radius (default 50 mm)
%        figure on: plots parameters, save: save a plot of parameter in
%        directory of motion parameters; off: nothing
%
% OUTPUT: FD framewise displacement
%         RMS the root mean square displacement 
%         
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

%% input
current = pwd;
radius  = 50;
fig     = 'save';
FD      = NaN;
RMS     = NaN;
motion  = NaN;

if nargin == 0
    [filename,filepath,sts] = uigetfile('*.txt','select realignement parameters');
    if sts == 0
        return
    else
        realignment_file = [filepath filesep filename];
    end
elseif nargin >=1
    realignment_file = varargin{1};
    [filepath,filename]=fileparts(realignment_file);
    for v=2:nargin
        if strcmpi(varargin{v},'Radius')
            radius = varargin{v+1};
        elseif strcmpi(varargin{v},'Figure')
            fig = varargin{v+1};
        end
    end
end

%% read
motion            = load(realignment_file); % load data
motion            = motion(:,1:6); % make sure to pick only the initial 6 if another file is used
motion            = spm_detrend(motion,1); % de-mean and detrend
motion(:,[4 5 6]) = motion(:,[4 5 6]).*radius; % angle in radian by average head size = displacement in mm

%% compute
D            = diff(motion,1,1); % 1st order derivative
D            = [zeros(1,6); D];  % set first row to 0
FD           = sum(abs(D),2); % framewise displacement a la Powers
RMS          = sqrt(mean(detrend(D).^2,2)); % root mean square for each column a la Van Dijk
FD_outliers  = spmup_comp_robust_outliers(FD);
RMS_outliers = spmup_comp_robust_outliers(RMS);

%% figure
if strcmpi(fig,'on') || strcmpi(fig,'save')
    figure('Name','Motion and displacement')
    set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1); plot(motion(:,[1 2 3]),'LineWidth',3); axis tight; box on; grid on; title('translation')
    xlabel('Volumes'); ylabel('displacement in mm')
    subplot(2,2,2); plot(motion(:,[4 5 6]),'LineWidth',3); axis tight; box on; grid on; title('rotation')
    xlabel('Volumes'); ylabel('displacement in mm')
    subplot(2,2,3); plot(FD,'LineWidth',2); axis tight; box on; grid on; title('framewise displacement')
    hold on; tmp = FD_outliers.*FD; tmp(tmp==0)=NaN; plot(tmp,'ro','LineWidth',3);
    xlabel('Volumes'); ylabel('absolute displacement in mm')
    subplot(2,2,4); plot(RMS,'LineWidth',2); axis tight; box on; grid on; title('root mean squares')
    hold on; tmp = RMS_outliers.*RMS; tmp(tmp==0)=NaN; plot(tmp,'ro','LineWidth',3);
    xlabel('Volumes'); ylabel('average displacement in mm')
    
    if strcmpi(fig,'save')        
        try
            print (gcf,'-dpdf', '-bestfit', fullfile(filepath,[filename(4:end) '_displacement.pdf']));
        catch
            print (gcf,'-dpdf', fullfile(filepath,[filename(4:end) '_displacement.pdf']));
        end
        close(gcf)
        cd(current)
    end
end

