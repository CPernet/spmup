function [r_course, outliers] = spmup_volumecorr(fmridata,skip)

% compute the Pearson correlation between volumes of the time series 
% and return possible outliers among correlated volumes (this assumes of
% course that r is not biased by outlying voxels between volumes)
%
% FORMAT [r_course, outliers] = spmup_volumecorr(fmridata,skip)
%
% INPUT fmridata is the name of the multiple 3D files or 4D file of fmri data
%       skip (optional) 
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

if nargin == 0
    [fmridata,sts]= spm_select(Inf,'image','Select fMRI data');
    if sts == 0; disp('selection aborded'); return; end
    if size(fmridata,1) == 1 && strcmp(fmridata(length(fmridata)-1:end),',1')
        fmridata = fmridata(1:length(fmridata)-2); % in case picked 4D put left ,1
    end
end

skipn = 0;
if exist('skip','var')
    skipn = skip; clear skip;
end

% check fmridata is 4D
if iscell(fmridata)
    for v=skipn+1:size(fmridata,1)
        Vfmri(v) =spm_vol(fmridata{v});
    end
else
    Vfmri = spm_vol(fmridata);
end

if sum(size(Vfmri)) == 2
    error('fMRI data must be time series')
end

% simply compute the Pearson correlation  
data = spm_read_vols(Vfmri);
data = reshape(data,[prod(Vfmri(1).dim),size(data,4)]);
r_course = corr(data); diag = 1:length(r_course)-1;
r_course = r_course(sub2ind(size(r_course),diag,(diag+1)));

% get the outliers from the r_course
if  nargout == 2
    % interquartile range
    y=sort(r_course);
    j=floor(length(r_course)/4 + 5/12);
    g=(length(r_course)/4)-j+(5/12);
    ql=(1-g).*y(j)+g.*y(j+1); % lower quartile
    k=length(r_course)-j+1;
    qu=(1-g).*y(k)+g.*y(k-1); % higher quartile
    value=qu-ql; % inter-quartile range
    
    % robust outliers
    M = median(r_course);
    k=(17.63*length(r_course)-23.64)/(7.74*length(r_course)-3.71); % Carling's k
    outliers=r_course<(M-k*value) | r_course>(M+k*value);
    
    % quick plot
%     figure; plot(r_course,'LineWidth',2); axis tight; grid on; title('Correlation time course','Fontsize',12)
%     xlabel('scans'); ylabel('r'); hold on; plot(r_course.*outliers,'or','LineWidth',2); 
%     axis([1 length(r_course(:)) min(r_course(:)) max(r_course(:))])

end

