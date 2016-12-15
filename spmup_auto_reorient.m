function RM = spmup_auto_reorient(p,which_image)

% FORMAT: reorientation_matrix = spmup_auto_reorient(p,which_image)
%
% IMPUT: p is a cell string of images to process (see spm_select)
%        which_image is a integer for wich image to obtain M
%
% OUTPUT: RM is the reorientation matrix obtain from affine reg to the
%         template - assuming the T1 and EPI are pretty much aligned, input
%         the T1 and use M into the realign utility buit in SPM to reorient
%         the EPI data (now all is set properly to reaign, segment, nornmalize)
%
% the function realigns each images with the template using affine transform
% this is the same idea as setting the AC point% manually except that
% (0,0,0) is now the same as the template (a priori better) -
%
% based on auto_reorient.m from Kiyotaka Nemoto
% http://www.nemotos.net/?p=17
%
% Cyril Pernet June 2015
% ------------------------
% spmup copyright

if ~nargin
    [p,sts] = spm_select(Inf,'image','Select the image(s) to re-orient');
    if ~sts, return; end
elseif nargin == 1
    which_image = 1;
end

p = cellstr(p);
if numel(p) == 1; which_image = 1; end
vg = spm_vol(fullfile(spm('Dir'),'canonical','avg152T1.nii'));
tmp = [tempname '.nii'];

fprintf('spmup: reorienting image %g \n',which_image);
spm_smooth(p{which_image},tmp,[12 12 12]);
vf = spm_vol(tmp);
M  = spm_affreg(vg,vf,struct('regtype','rigid'));
[u,s,v] = svd(M(1:3,1:3));
M(1:3,1:3) = u*v';
RM = M;
N  = nifti(p{which_image});
N.mat = RM*N.mat;
create(N);
spm_unlink(tmp);

for i=1:numel(p)
    if i ~= which_image
        fprintf('spmup: applying to image %g \n',i);
        N  = nifti(p{i});
        N.mat = RM*N.mat;
        create(N);
    end
end
spm_unlink(tmp);

