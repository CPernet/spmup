% ut_spmup_anatQA

% reproduce the code but uses simple 2D image ensuring the output follows
% what is expected

%% Case 1 - good qality image 
% the image is 64*64 with the brain a 32*32 box in the middle (itself a
% 16*16GM and WM)

myimage = zeros(64,64); % image
mybrain = abs(randn(32,32)*100);
myimage(17:48,17:48) = mybrain; % gray matter
myimage(25:40,25:40) = myimage(25:40,25:40)+50; %white matter

% brain_mask = (smooth3(spm_read_vols(GrayV),'box',25)+smooth3(spm_read_vols(WhiteV),'box',25))>0;
% [x,y,z] = ind2sub(AnatV.dim,find(brain_mask==0));
brain_mask = find(myimage==0);
data = sort(myimage(brain_mask)); 
% std_nonbrain = std(data(data<(median(data)+iqr(data)/2))); % in case we have high values (ie brain, ghost) remove top 25%
std_nonbrain = 1;

% [x,y,z] = ind2sub(AnatV.dim,find(spm_read_vols(GrayV) > 0.6));
% dataGM = spm_get_data(AnatV,[x y z]');
% meanGM = mean(data);
GM1 = myimage(17:24,17:24); GM2 = myimage(40:48,40:48);
dataGM = [GM1(:); GM2(:)];
meanGM = mean(dataGM);

% [x,y,z] = ind2sub(AnatV.dim,find(spm_read_vols(WhiteV) > 0.6));
% dataWM = spm_get_data(AnatV,[x y z]');
dataWM = myimage(25:40,25:40);
dataWM = dataWM(:); meanWM = mean(dataWM);

SNR  = ((meanGM+meanWM)/2)/std_nonbrain;
CNR  = (meanWM-meanGM)/std_nonbrain;
FBER = var([dataGM;dataWM])/std_nonbrain^2;
Bmax = sqrt(sum(myimage(:).^2));
EFC = nansum((myimage(:)./Bmax).*abs(log((myimage(:)./Bmax))));

figure; subplot(1,3,1);
imagesc(myimage); title(sprintf('Clean image \n SNR=%g CNR=%g \n FBER=%g EFC=%g',SNR,CNR,FBER,EFC));

%% Case 2 - add general white noise 
% the image is 64*64 with the brain a 32*32 box in the middle (itself a
% 16*16GM and WM)

myimage = randn(64,64)*10;
myimage(17:48,17:48) = mybrain; % gray matter
myimage(25:40,25:40) = myimage(25:40,25:40)+50; %white matter
data = sort(myimage(brain_mask)); 
std_nonbrain = std(data(data<(median(data)+iqr(data)/2))); % in case we have high values (ie brain, ghost) remove top 25%
GM1 = myimage(17:24,17:24); GM2 = myimage(40:48,40:48);
dataGM = [GM1(:); GM2(:)];
meanGM = mean(dataGM);
dataWM = myimage(25:40,25:40);
dataWM = dataWM(:); meanWM = mean(dataWM);
SNR  = ((meanGM+meanWM)/2)/std_nonbrain;
CNR  = (meanWM-meanGM)/std_nonbrain;
FBER = var([dataGM;dataWM])/std_nonbrain^2;
Bmax = sqrt(sum(myimage(:).^2));
EFC = nansum((myimage(:)./Bmax).*abs(log((myimage(:)./Bmax))));

subplot(1,3,2);
imagesc(myimage); title(sprintf('Adding White noise (10) \n SNR=%g CNR=%g \n FBER=%g EFC=%g',SNR,CNR,FBER,EFC));

%% Case 3 - add motion, i.e. blurring 
% the image is 64*64 with the brain a 32*32 box in the middle (itself a
% 16*16GM and WM)

mybrain = myimage(17:48,17:48);

myimage = myimage + noise;
data = sort(myimage(brain_mask)); 
std_nonbrain = std(data(data<(median(data)+iqr(data)/2))); % in case we have high values (ie brain, ghost) remove top 25%
GM1 = myimage(17:24,17:24); GM2 = myimage(40:48,40:48);
dataGM = [GM1(:); GM2(:)];
meanGM = mean(dataGM);
dataWM = myimage(25:40,25:40);
dataWM = dataWM(:); meanWM = mean(dataWM);
SNR  = ((meanGM+meanWM)/2)/std_nonbrain;
CNR  = (meanWM-meanGM)/std_nonbrain;
FBER = var([dataGM;dataWM])/std_nonbrain^2;
Bmax = sqrt(sum(myimage(:).^2));
EFC = nansum((myimage(:)./Bmax).*abs(log((myimage(:)./Bmax))));

subplot(1,3,2);
imagesc(myimage); title(sprintf('Adding White noise (10) \n SNR=%g CNR=%g \n FBER=%g EFC=%g',SNR,CNR,FBER,EFC));
