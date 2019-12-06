function spmup_spectral

% fMRI spectral analysis - slice by slice
% Cyril Pernet - University of Edinburgh 10 February 2017

%%   get the data

clear var; warning off; current = pwd;
path = uigetdir(pwd,'select your working directory');
cd(path); defaults = spm_get_defaults; global defaults;
[P, sts]= spm_select(Inf,'image','Select Image time series');
if sts == 0 || isempty(P) == 1
    disp('no item selected')
    return
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

%% computye the image spectral power slice by slice
f=zeros(nbimage,zmax,256/2+1);
pow = zeros(256,256,zmax-1);
for j=1:nbimage
    for i=1:zmax-1
        im = squeeze(Image(:,:,i,j)); im(isnan(im)) = 0;
        % figure; subplot(1,3,1); imagesc(im)
        % make the image square either by zero padding or truncate
        N = min(size(im));
        index = (max(size(im)) - N) / 2;
        im = im((1+index):size(im,1)-index,:);
        imf2=fftshift(fft2(im,256,256)); % zero-pads im to be 256-by-256
        freq =-256/2:256/2-1; subplot(2,2,2);
        impf=abs(imf2).^2; pow(:,:,i) = impf;
        % subplot(1,3,2); imagesc(freq,freq,log10(impf));
        [X Y]=meshgrid(freq,freq); % get coordinates of the power spectrum image
        [theta rho]=cart2pol(X,Y); % equivalent in polar coordinates
        rho=round(rho);
        for r=0:256/2
            ii{r+1}=find(rho==r); % for each freq return the location in the polar
            % array ('find') of the corresponding frequency
            f(j,i,r+1)=mean(impf(ii{r+1})); % average power values of all the same freq
        end
        % subplot(1,3,3); freq2=0:256/2; loglog(freq2,squeeze(f(j,i,:))); title('frequency spectrum','Fontsize',14); axis tight
        
%         figure % fig for odd vs even slices
%         D = mean(pow(:,:,1:2:zmax-1),3)-mean(pow(:,:,2:2:zmax-1),3);
%         imagesc(freq,freq,log10(abs(D))); draw now
    end
end

%% figure
avgf = squeeze(mean(f,1)); % avg over time series
figure; freq2=0:256/2; loglog(freq2,avgf([1:2:zmax],:),'k');
hold on; loglog(freq2,avgf([2:2:zmax],:),'r'); axis tight

