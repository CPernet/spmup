function [rP,rC]=spmup_corr(image1,image2,mask,figout,threshold)

% Computes the Pearson (vector-wise) and concordance correlations
% between image1 and image2 for data included in the mask
%
% FORMAT = [rP,rC]=spmup_corr(image1,image1,mask,threshold)
%
% INPUT image1 is the filename of an image (see spm_select)
%       image2 is the filename of an image (see spm_select)
%       mask   is the filename of a mask in same space as image1 and image2
%       figout 1/0 (default) to get correlation figure out
%       threshold (optional) if mask is not binary, threshold to apply
%
% OUTPUT rP is the Pearson correlation coefficient
%        rC is the concordance correlation coefficient
%
% Concordance correlation is more useful for reliability because it estimates
% how much variation from the 45 degree line we have (by using the cov)
% see Lin, L.I. 1989. A Corcordance Correlation Coefficient to Evaluate
% Reproducibility. Biometrics 45, 255-268.
%
% Cyril Pernet
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2017

%% check inputs
spm('defaults', 'FMRI');
V1 = spm_vol(image1);
V2 = spm_vol(image2);
if any(V1.dim ~= V2.dim)
    error('input images are of different dimensions')
end

M = spm_vol(mask);
if any(M.dim ~= V1.dim)
    error('mask and input image are of different dimensions')
end
mask = spm_read_vols(M);
if exist('threshold','var')
    mask = mask > threshold;
end

%% Pearson correlation
[x,y,z] = ind2sub(M.dim,find(mask));
X = [spm_get_data(V1,[x y z]'); spm_get_data(V2,[x y z]')]';
X(isnan(X(:,1)),:) = []; % clean up if needed
rP = sum(detrend(X(:,1),'constant').*detrend(X(:,2),'constant')) ./ ...
    (sum(detrend(X(:,1),'constant').^2).*sum(detrend(X(:,2),'constant').^2)).^(1/2);

if figout == 1
    figure; scatter(X(:,1),X(:,2),50); grid on
    xlabel('img1','FontSize',14); ylabel('img2','FontSize',14);
    title(['corr =' num2str(rP)],'FontSize',16); box on;set(gca,'Fontsize',14)
    h=lsline; set(h,'Color','r','LineWidth',4);
end


%% Concordance
if nargout == 2
    S = cov(X,1); Var1 = S(1,1); Var2 = S(2,2); S = S(1,2);
    rC = (2.*S) ./ (Var1+Var2+((mean(X(:,1)-mean(X(:,2)).^2))));
end


