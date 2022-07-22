% ut_spmup_anatQA
%
% call spmup_anatQA testing with simulated data
% that means all the data in GM and WM masks are used corresponding to the
% values generated therein 
%
% Cyril Pernet 
% -------------------------------------------------------------------------
% Copyright (C) spmup team 

%% load up the SPM tpm

clearvars
spmroot = fileparts(which('spm'));
for n=6:-1:1
    P(n,:) = [spmroot filesep 'tpm' filesep 'TPM.nii,' num2str(n)];
end
V = spm_vol(P);

I_brain = spm_read_vols(V);
I_brain = sum(I_brain(:,:,:,1:5),4); % add all tissue types but last

%% make the image with known signal

% background
I_Bckgd         = spm_read_vols(V(6))>0.1; % threshol prob at 0.1
index           = find(I_Bckgd);
noise           = randn(length(index),1); % generate noise for background
I_brain(index)  = noise;
clear I_Bckgd

% gray matter and white matter
tmp_GM          = spm_read_vols(V(1)); 
tmp_WM          = spm_read_vols(V(2)); 
I_GM            = tmp_GM > tmp_WM;
I_WM            = tmp_WM > tmp_GM; % make gray/white mutually exclusive
clear tmp_GM tmp_WM 

index1          = find(I_GM);
GM              = randn(length(index1),1)*50; % generate data for GM
I_brain(index1) = GM;

index2          = find(I_WM);
WM              = randn(length(index2),1)*300; % generate data for WM
I_brain(index2) = WM;

tmpdir = fullfile(fileparts(fileparts(which('spmup_anatQA.m'))),['unit_test' filesep 'tmp_anatQA']); 
mkdir(tmpdir); W=V(1); W.fname = [tmpdir filesep 'brain.nii'];
W.descrip = 'fake brain for QA'; spm_write_vol(W,I_brain);

%% Compute QA from generated data
%        - SNR : the signal-to-Noise Ratio, ie the mean intensity within gray and white matter divided
%                by the standard deviation of the values outside the brain. Higher values are better.
%        - CNR : the Contrast to Noise Ratio, i.e. the mean of the white matter intensity values minus the mean 
%                of the gray matter intensity values divided by the standard deviation of the values outside the
%                brain. Higher values are better.
%        - FBER: Foreground to Background Energy Ratio, i.e. the variance of voxels in grey and white matter 
%                divided by the variance of voxels outside the brain. Higher values are better.
%        - EFC : Entropy Focus Criterion, i.e. the entropy of voxel intensities proportional to the maximum 
%                possibly entropy for a similarly sized image. Indicates ghosting and head motion-induced blurring. 
%                Lower values are better. See <http://ieeexplore.ieee.org/document/650886/>

% see what values we get using spmup_anatQA
[anatQA,int_data] = spmup_anatQA(W.fname,P(1,:),P(2,:));

% compute from simulated data, except noise because it's compute differently 
SNR  = ((mean(GM)+mean(WM))/2) / int_data.std_nonbrain;
CNR  = (mean(WM)-mean(GM)) / int_data.std_nonbrain;
FBER = var([GM ;WM]) / int_data.std_nonbrain^2;
Bmax = sqrt(sum(I_brain(:).^2));
EFC  = nansum((I_brain(:)./Bmax).*abs(log((I_brain(:)./Bmax))));

if ~any([single(anatQA.SNR)==single(SNR),...
        single(anatQA.CNR)==single(CNR),...
        single(anatQA.FBER)==single(FBER),...
        single(anatQA.EFC)==single(EFC)])
    warndlg('simulated data and measures obtained are not the same !!!')
else
    disp('all good SNR, CNR, FBER and EFC values out are equal to simulated ones')
end

% still check noise levels are comparable
up         = median(noise)+iqr(noise); 
down       = median(noise)-iqr(noise);
index      = logical((noise<up).*(noise>down)); 
std_noise1 = std(noise(index)); 
std_noise2 = std(noise); 
if int_data.std_nonbrain>std_noise1 && int_data.std_nonbrain<std_noise2
    disp('noise_level estimate ok as well')
else
    warndlg('simulated data and noise estimate out of bounds')
end

% figure; imagesc(squeeze(I_brain(:,:,round(V(1).dim(3)/2)))); colormap('gray')
% title(sprintf('SNR=%g CNR=%g \n FBER=%g EFC=%g',SNR,CNR,FBER,EFC));
rmdir(tmpdir,'s')






