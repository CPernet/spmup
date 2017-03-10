function [volume_outliers, slice_outliers] = spmup_spatialcorr_qa(P)

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

clear all; warning off; current = pwd;
path = uigetdir(pwd,'select your working directory');
cd(path); defaults = spm_get_defaults; global defaults; 
if nargin == 0
    [P, sts]= spm_select(Inf,'image','Select Image time series');
    if sts == 0 || isempty(P) == 1
        disp('no item selected')
        return
    end
end

V = spm_vol(P);
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
nbimage = size(V,1);


%% ---------------
% get a mean and std images
% compute a mask to 'see' the edges

disp('computing mean, std, mask')
mean_image = sum(Image,4) / nbimage ;
std_image = std(Image,0,4) ;

% compute a mask
% assume threshold as the 2st bin of the histogramm

for n=1:nbimage
    binary_img(:,:,:,n)=Image(:,:,:,n);
    for z = 1:zmax
        [N,X]=hist(binary_img(:,:,z));
        binary_img(:,:,z,n)=Image(:,:,z,n) > X(2);
    end
end

mask = sum(binary_img,4);
mask = (mask == nbimage);

% save info
fprintf('writing mean, std, and mask images to disk \n output in %s\n',path)
mkdir('fMRI_info'); cd fMRI_info
path = pwd; Info_img = V(1);
name = '/mean_image.nii';
Info_img.fname = sprintf('%s%s',path,name);
Info_img.description = 'mean_image';
spm_write_vol(Info_img,mean_image);
name = '/std_image.nii';
Info_img.fname = sprintf('%s%s',path,name);
Info_img.description = 'std_image';
spm_write_vol(Info_img,std_image);
name = '/mask.nii';
Info_img.fname = sprintf('%s%s',path,name);
Info_img.description = 'mask image';
spm_write_vol(Info_img,mask);

% --------------------------------------
% display

figure
for i=1:zmax
    subplot(1,2,1)
    image(rot90(mean_image(:,:,i).*mask(:,:,i)))
    title('mean')
    subplot(1,2,2)
    image(rot90(std_image(:,:,i).*mask(:,:,i)))
    title('std')
    colormap('jet')
    drawnow
    pause(0.2)
end


%% ---------------
% check the different slices / volume per volumes
% simply use a linear correlation slice per slice
% a r distribution is then obtained per location for
% the time series and one checks for outlier

for j=1:nbimage
    for i=1:zmax-1
        if mean(mean(mask(:,:,i+1))) ~= 0
            r(i,j) = corr2(Image(:,:,i,j),Image(:,:,i+1,j));
        else
            r(i,j) = NaN;
        end
    end
end

figure; imagesc(r); ylabel('slices'); xlabel('volume nb')
title('Correlation between slices per volume','FontSize',14)
drawnow


data = nanmean(r,1); y=sort(data);
j=floor(length(data)/4 + 5/12);
g=(length(data)/4)-j+(5/12);
ql=(1-g).*y(j)+g.*y(j+1); % lower quartile
k=length(data)-j+1;
qu=(1-g).*y(k)+g.*y(k-1); % higher quartile
value=qu-ql; % inter-quartile range
M = median(data);
k=(17.63*n-23.64)/(7.74*n-3.71); % Carling's k
volume_outliers=data<(M-k*value) | data>(M+k*value);
if sum(volume_outliers ~= 0)
    fprintf('volume %g is potiential outlier\n',find(volume_outliers));
else
    disp('no volume outliers detected');
end




%% average neighbouring slices
outliers = find(volume_outliers);
if sum(diff(outliers) == 1) == 0
    answ = input('do you want to ''repair'' bad slices using averaging? [y/n] ','s');
    if strcmpi(answ,'y')
        for im=1:length(outliers)
            bad_slices = r(:,outliers(im))<(M-k*value);
            Image(:,:,bad_slices,outliers(im)) = (Image(:,:,bad_slices,outliers(im)-1) + ...
                Image(:,:,bad_slices,outliers(im)+1))/2;
        end
        
        for im=1:size(V,1) 
            [p,f,e]=fileparts(V(im).fname);
            V(im).fname = [p filesep 'repaired_' f e];
            spm_write_vol(V(im),Image(:,:,:,im));
        end
    end
else
    disp('several consecutive volumes are seen as outliers, sorry cannot repair these data')
end
    
%% check outliers between volumes per slice
for j=1:nbimage-1
    for i=1:zmax
        if mean(mean(mask(:,:,i))) ~= 0
            r(i,j) = corr2(Image(:,:,i,j),Image(:,:,i,j+1));
        else
            r(i,j) = NaN;
        end
    end
end

figure; imagesc(r)
title('Correlation between volumes per slice ','FontSize',14)
ylabel('slices'); xlabel('volume nb')
drawnow

data = nanmean(r,2); y=sort(data);
j=floor(length(data)/4 + 5/12);
g=(length(data)/4)-j+(5/12);
ql=(1-g).*y(j)+g.*y(j+1); % lower quartile
k=length(data)-j+1;
qu=(1-g).*y(k)+g.*y(k-1); % higher quartile
value=qu-ql; % inter-quartile range
M = median(data);
k=(17.63*n-23.64)/(7.74*n-3.71); % Carling's k
slice_outliers=data<(M-k*value) | data>(M+k*value);
if sum(slice_outliers ~= 0)
    fprintf('slice %g is potiential outlier\n',find(slice_outliers));
else
    disp('no slice outliers detected');
end





