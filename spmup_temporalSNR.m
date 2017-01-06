function tSNR = spmup_temporalSNR(time_series,masks,figout)

% Computes the temporal SNR of the time_series input in the different
% compartments provided by the masks images = mean signal / std over time
% The routine recapitulates tSNR as described in Thomas Liu  (2016) 
% Noise contributions to the fMRI signal: An overview NeuroImage, 343, 
% 141-151 <http://dx.doi.org/10.1016/j.neuroimage.2016.09.008>
% 
% FORMAT tSNR = spmup_temporalSNR(time_series,masks,plot)
%
% INPUT time_series: a cell array of file names (see spm_select)
%       masks: a cell array of (non binary) tissue class images (GM, WM, 
%       CSF in that order) in the same space as the time series (i.e. 
%       typically from the anatomical coregistered to the mean EPI)
%       plot 0/1 if you want to have a plot of tSNR per ROI size
%
% OUTPUT tSNR is a structure with the following fields:
%            .GM: mean GM signal / std over time (estimate BOLD from GM>(WM+CSF))
%            .WM:  mean WM signal / std over time (estimate non-BOLD from WM>(GM+CSF))
%            .CSF: mean CSF signal / std over time (estimate non-BOLD from CSF>(GM+WM))
%            .Background:  mean signal outside mask (GM+WM+CSF) / std over time
%            .average (tSNR): mean signal / sqrt(std(GM)+std(WM+CSF)+std(Background))
%            .roi: tSNR for increased ROI (from masks 0.95 to 0.1) ~linear function of srqrt(nb voxels)
%            .image (SNR0):  mean signal inside mask / std outside mask over time
%            .signal_mean: [std(GM)+std(WM+CSF)] / sqrt((SNR0^2/tSNR- 1)/SNR0^2)
%            Since tSNR = SNR0^2 / (1+L^2*SNR0^2), we have L^2 = (SNR0^2 /tSNR - 1) / SNR0^2
%            and std(GM)+std(WM+CSF) = L*Smean 
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

%% check inputs
if nargin == 2
    figout = 1;
end

V = spm_vol(time_series);
if size(V,1) < 10
    error('there is less than 10 images in your time series ??')
end

VM = spm_vol(masks);
if size(VM,1) ~= 3
    error(['3 masks files expected, ' num2str(size(VM,1)) ' detected - check input file'])
end

    
%% Compute relative masks
GM = spm_read_vols(VM(1));
WM = spm_read_vols(VM(2));
CSF = spm_read_vols(VM(3));

brain_mask = (smooth3(GM)+smooth3(WM)+smooth3(CSF))>0;
% figure; for z=1:size(brain_mask,3); imagesc(squeeze(brain_mask(:,:,z))); pause; end
GM = GM.*(GM>0.2); WM = WM.*(WM>0.2); CSF = CSF.*(CSF>0.2); % baseline prob 20%
GM = GM.*(GM>(WM+CSF)); % figure; rst_hist(GM(:))
WM = WM.*(WM>(GM+CSF)); % figure; rst_hist(WM(:))
CSF = CSF.*(CSF>(WM+GM)); % figure; rst_hist(CSF(:))

%% in masks tSNR
[x,y,z] = ind2sub(size(GM),find(GM));
data = spm_get_data(V,[x y z]');
stdGM = mean(std(data,1));
tSNR.GM =  mean(mean(data,1)) /stdGM;

[x,y,z] = ind2sub(size(WM),find(WM));
data = spm_get_data(V,[x y z]');
stdWM = mean(std(data,1));
tSNR.WM =  mean(mean(data,1)) /stdWM;

[x,y,z] = ind2sub(size(CSF),find(CSF));
data = spm_get_data(V,[x y z]');
stdCSF = mean(std(data,1));
tSNR.CSF =  mean(mean(data,1)) /stdCSF;

%% Background
[x,y,z] = ind2sub(size(brain_mask),find(brain_mask ~= 1));
data = spm_get_data(V,[x y z]');
stdBackground = mean(std(data,1));
tSNR.Background =  mean(mean(data,1)) /stdBackground;

%% average (tSNR)
[x,y,z] = ind2sub(size(WM),find(WM+CSF));
data = spm_get_data(V,[x y z]');
stdWMCSF = mean(std(data,1)); % presumably non BOLD

[x,y,z] = ind2sub(size(brain_mask),find(GM+WM+CSF+(brain_mask ~= 1)));
data = spm_get_data(V,[x y z]'); % the whole image or so
tSNR.average = mean(mean(data,1)) / sqrt(stdGM^2+stdWMCSF^2+stdBackground^2);

%% SNR0
[x,y,z] = ind2sub(size(brain_mask),find(brain_mask));
data = spm_get_data(V,[x y z]');
tSNR.image = mean(mean(data,1)) / sqrt(stdGM^2+stdWMCSF^2);

%% signal
L2 = (tSNR.image^2 /tSNR.average - 1) / tSNR.image^2;
tSNR.signal_mean = sqrt(stdGM^2+stdWMCSF^2) / sqrt(L2);


%% per ROI (absolute masking)
GM = spm_read_vols(VM(1));
WM = spm_read_vols(VM(2));
CSF = spm_read_vols(VM(3));

index = 1;
for p=0.95:-0.05:0.1
  
    [x,y,z] = ind2sub(size(GM),find(GM>=p));
    data = spm_get_data(V,[x y z]');
    stdGM = mean(std(data,1)); 

   [x,y,z] = ind2sub(size(WM),find(WM>p+CSF>p));
   data = spm_get_data(V,[x y z]');
   stdWMCSF = mean(std(data,1)); 

   ROI = (GM>p+WM>p+CSF>p);
   % figure; for z=1:size(ROI,3); imagesc(squeeze(ROI(:,:,z))); pause; end
   [x,y,z] = ind2sub(size(brain_mask),find(ROI+(brain_mask~=1)));
   data = spm_get_data(V,[x y z]'); % the whole image or so
   tSNR.roi.value(index) = mean(mean(data,1)) / sqrt(stdGM^2+stdWMCSF^2+stdBackground^2);
   tSNR.roi.size(index) = sum(ROI(:));
   index = index + 1;
end
B = pinv([sqrt(tSNR.roi.size)' ones(18,1)])*tSNR.roi.value';
model = [sqrt(tSNR.roi.size)' ones(18,1)]*B;
tSNR.roi.slope = B(1);

if figout ~=0
    figure('Name','spmup tSNR');
    plot(sqrt(tSNR.roi.size),[sqrt(tSNR.roi.size)' ones(18,1)]*B,'LineWidth',3);
    hold on; plot(sqrt(tSNR.roi.size),tSNR.roi.value,'ro','LineWidth',2);
    axis tight; box on; grid minor; ylabel('temporal SNR','FontSize',12)
    xlabel('sqrt of the number of in brain voxels used','FontSize',12)
    title(['linear fit temporal SNR/in brain ROI size RMSE=' num2str(sqrt(mean(model - tSNR.roi.value')))],'FontSize',12)
end

    