function [volume_outliers, slice_outliers] = spmup_spatialcorr_qa(P,spath)

% look at the time series per slice, and compute correlations
% this allows testing that images are alike and spots drop out easily
% also create the mean image, sd, and binary mask 
%
% FORMAT: [volume_outliers, slice_outliers] = fMRI_info(P)
% INPUT:  P is a cell array of image names
% OUTPUT: volume_outliers indexes volumes with an average correlation between
%         slices bigger or lower than the other volumes
%         slice_outliers indexes slices with an average correlation between
%         volumes bigger or lower than the other slices
%
% The outlier detection is the adjusted Carling's box-plot rule
% ie outside the bound of median+/- k*IQR, 
% k = k=(17.63*n-23.64)/(7.74*n-3.71)
% see Carling, K. (2000). Resistant outlier rules and the non-Gaussian case. 
% Stat. Data Anal. 33, 249–258. 
% <http://www.sciencedirect.com/science/article/pii/S0167947399000572>
% The Inter-Quartile-Range is computed the ideal fourths 
% see D.C. Hoaglin, B. Iglewicz (1987) Fine-tuning some resistant rules for
% outlier labelling. J. Amer. Statist. Assoc., 82 , 1147–1149
% <http://www.tandfonline.com/doi/abs/10.1080/01621459.1986.10478363>
%
% Relies on SPM functions (see http://www.fil.ion.ucl.ac.uk/spm/)
% Cyril Pernet v4 - 25 Feb 2016, The University of Edinburgh
% ---------------------------------------------------------------

%%   get the data

warning off; current = pwd;
defaults = spm_get_defaults; global defaults; 
if nargin == 0
    [P, sts]= spm_select(Inf,'image','Select Image time series');
    if sts == 0 || isempty(P) == 1
        disp('no item selected')
        return
    end
    spath = uigetdir(pwd,'select your working directory');
end

cd(spath);
if iscell(P)
    for v=1:size(P,1)
        V(v) =spm_vol(P{v});
    end
else
    V = spm_vol(P);
end
N = numel(V);
Image = zeros([V(1).dim(1:3),N]);
for i=1:N
    for p=1:V(1).dim(3)
        Image(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
    end
end
xmax  = V(1).dim(1);
ymax  = V(1).dim(2);
zmax  = V(1).dim(3);


%% ---------------
% check the different slices / volume per volumes
% simply use a linear correlation slice per slice
% a r distribution is then obtained per location for
% the time series and one checks for outlier

for j=1:N
    for i=1:zmax-1
        if mean(mean(mask(:,:,i+1))) ~= 0
            r(i,j) = corr2(Image(:,:,i,j),Image(:,:,i+1,j));
        else
            r(i,j) = NaN;
        end
    end
end


figure; 
subplot(2,2,1); imagesc(r); ylabel('slices'); xlabel('volume nb')
title('Correlation between slices per volume','FontSize',14)

r = nanmean(r,1); y=sort(r);
j=floor(length(r)/4 + 5/12);
g=(length(r)/4)-j+(5/12);
ql=(1-g).*y(j)+g.*y(j+1); % lower quartile
k=length(r)-j+1;
qu=(1-g).*y(k)+g.*y(k-1); % higher quartile
value=qu-ql; % inter-quartile range
M = median(r);
k=(17.63*length(r)-23.64)/(7.74*length(r)-3.71); % Carling's k
volume_outliers=r<(M-k*value) | r>(M+k*value);
if sum(volume_outliers ~= 0)
    fprintf('volume %g is potiential outlier\n',find(volume_outliers));
else
    disp('no volume outliers detected');
end

subplot(2,2,3); plot(r,'LIneWidth',3); ylabel('average r'); xlabel('volume nb')
box on; grid on;

%% check outliers between volumes per slice
for j=1:N
    for i=1:zmax
        if mean(mean(mask(:,:,i))) ~= 0
            r(i,j) = corr2(Image(:,:,i,j),Image(:,:,i,j+1));
        else
            r(i,j) = NaN;
        end
    end
end

subplot(2,2,2); imagesc(r)
title('Correlation between volumes per slice ','FontSize',14)
ylabel('slices'); xlabel('volume nb')

r = nanmean(r,2); y=sort(r);
j=floor(length(r)/4 + 5/12);
g=(length(r)/4)-j+(5/12);
ql=(1-g).*y(j)+g.*y(j+1); % lower quartile
k=length(r)-j+1;
qu=(1-g).*y(k)+g.*y(k-1); % higher quartile
value=qu-ql; % inter-quartile range
M = median(r);
k=(17.63*length(r)-23.64)/(7.74*length(r)-3.71); % Carling's k
slice_outliers=r<(M-k*value) | r>(M+k*value);
if sum(slice_outliers ~= 0)
    fprintf('slice %g is potiential outlier\n',find(slice_outliers));
else
    disp('no slice outliers detected');
end

subplot(2,2,4); plot(r,'LIneWidth',3); ylabel('average r'); xlabel('volume nb')
box on; grid on;




