function [bad_subject] = find_group_outliers(SPM_file)

% routine that checks for outliers among subjects
% using a modified boxplot rule, each voxel is analyzed
% looking for outlier subjects and the proportion of outlying
% voxels is obtained per subjet allowing to decide, using again
% a modified box-plot rule, if a given subject has an abnormal
% amount of outlying voxels
%
% Cyril Pernet 15-07-2013

if nargin == 0
    [SPM_file,sts]=spm_select(1,'mat','select SPM.mat file');
    if sts == 0
        return
    end
    current = pwd;
end

% get the subject files
load(SPM_file);
V = spm_vol(SPM.xY.P);
for n=1:SPM.nscan
    Images(:,:,:,n) = spm_read_vols(V{n});
end

% get mask
[path,file,ext] = fileparts(SPM_file);
cd(path);
clear V;
Mask = spm_read_vols(spm_vol([pwd filesep 'mask.img']));

% apply mask on data
for n=1:SPM.nscan
    Images(:,:,:,n) = squeeze(Images(:,:,:,n)) .* Mask;
end

% now loop per voxel
Outliers = NaN(size(Images));
index = find(Mask);
disp('computing outliers per voxel')
for v=1:length(index)
    [x,y,z]=ind2sub(size(Mask),index(v));
    Outliers(x,y,z,:) = iqr_method(squeeze(Images(x,y,z,:)));
end

% finally check if one subject has an abnormal amount of outlying voxels
% indeed given the total number of voxel we will always find some outliers
% but it should be distributed equally across subjects
for n=1:SPM.nscan
    tmp = squeeze(Outliers(:,:,:,n));
    percent(n) = nansum(tmp(:)) / length(index);
end
bad_subject = find(iqr_method(percent,2));

% output
figure; hist(percent.*100); 
grid on; xlabel('percentage of outliying voxels');
ylabel('frequency'); title('histogram of the percentage of outlying voxels');

if isempty(bad_subject)
    disp('no outliers found')
else
    for n=1:length(bad_subject)
        [path,file,ext] = fileparts(cell2mat(SPM.xY.P(bad_subject(n))));
        fprintf('subject %s \n %s \n is likely an outlier',path,file)
    end
end

cd(current)

%% outlier detection
% -----------------------
function [I] = iqr_method(a,out)

% Returns a logical vector that flags outliers as 1s based 
% on the IQR methods described in Wilcox 2012 p 96-97.
%
% uses Carling's modification of the boxplot rule.
% An observation Xi is declared an outlier if:
% Xi<M-k(q2-q1) or Xi>M+k(q2-q1),
% where M is the sample median, and
% k=(17.63n-23.64)/(7.74n-3.71),
% where n is the sample size.
%
% INPUT a is a vector of data
%       out determines the thresholding 1 bilateral 2 unilateral
%
% OUTPUTS: I =     logical vector with 1s for outliers
%          value = IQR, the inter-quartile range 
%
% Cyril Pernet - adapted from Corr_toolbox
% ----------------------------------------

if nargin == 1
    out = 1;
end

a=a(:);n=length(a);
% inter-quartile range
j=floor(length(a)/4 + 5/12);
y=sort(a);
g=(length(a)/4)-j+(5/12);
q1=(1-g).*y(j)+g.*y(j+1);
k=length(a)-j+1;
q2=(1-g).*y(k)+g.*y(k-1);
value=q2-q1;
% outliers
M = median(a);
k=(17.63*n-23.64)/(7.74*n-3.71);
if out == 1
    I=a<(M-k*value) | a>(M+k*value);
else
    I= a>(M+k*value); % only reject data with a too high value
end
I = I+isnan(a);
