function RM = spmup_auto_reorient(P,which_image)

% FORMAT: reorientation_matrix = spmup_auto_reorient(p,which_image)
%
% IMPUT: p is a cell string of images to process (see spm_select)
%        which_image is a integer for wich image to obtain M
%
% OUTPUT: RM is the reorientation matrix obtain from affine reg to the
%         template - assuming the T1 and EPI are pretty much aligned, it
%         takes the T1 and use M into the realign utility buit in SPM to reorient
%         the EPI data (now all is set properly to segment and normalize)
%
% the function realigns each images with the template using affine transform
% this is the same idea as setting the AC point manually except that
% (0,0,0) is now the same as the template (a priori better) -
%
% based on auto_reorient.m from Kiyotaka Nemoto
% http://www.nemotos.net/?p=17
%
% Cyril Pernet June 2015
% ------------------------
% spmup copyright

if ~nargin
    [P,sts] = spm_select(Inf,'image','Select the image(s) to re-orient');
    if ~sts, return; end
elseif nargin == 1
    which_image = 1;
end

addpath([spm_file(which('spm'),'path') filesep 'toolbox' filesep 'OldNorm']);
V = spm_vol(P);    
if numel(V) == 1; which_image = 1; end
vg = spm_vol(fullfile(spm('Dir'),'canonical','avg152T1.nii'));
tmp = [tempname '.nii'];

fprintf('spmup: reorienting image %g \n',which_image);
if size(P,1) == 1 && numel(V) == 1 % 1 image only
    spm_smooth(P,tmp,[12 12 12]);
elseif numel(V) >1 % series
    try
        if numel(V{which_image}) == 1
            spm_smooth(P{which_image},tmp,[12 12 12]);
        else % numel(V) > 1 % multiple 3D
            spm_smooth([P ',' num2str(which_image)],tmp,[12 12 12]);
        end
    catch
        if numel(V) > 1 % multiple 3D
            spm_smooth([P ',' num2str(which_image)],tmp,[12 12 12]);
        else 
            spm_smooth(P(which_image),tmp,[12 12 12]);
        end
    end
end
vf = spm_vol(tmp);
M  = spm_affreg(vg,vf,struct('regtype','rigid'));
[u,~,v] = svd(M(1:3,1:3));
M(1:3,1:3) = u*v';
RM = M;
if size(P,1) == 1 && numel(V) == 1 % 1 image only
    N  = nifti(P);
elseif numel(V) >1 % series
    try
        if numel(V{which_image}) == 1
        N  = nifti(P{which_image});
    else % numel(V) > 1 % multiple 3D 
        N  = nifti([P ',' num2str(which_image)]);
        end
    catch
        if numel(V) > 1 % multiple 3D
            N  = nifti([P ',' num2str(which_image)]);
        else 
            N  = nifti(P(which_image));
        end
    end
end

N.mat = RM*N.mat;
create(N);

if size(P,1) >1
    for i=1:numel(V)
        if i ~= which_image
            fprintf('spmup: applying to image %g \n',i);
            N  = nifti(P{i});
            N.mat = RM*N.mat;
            create(N);
        end
    end
else
    if numel(V) > 1
        fprintf('spmup: applying to all other images \n');
    end
    N  = nifti(P);
    N.mat = RM*N.mat;
    create(N);
end
spm_unlink(tmp);


