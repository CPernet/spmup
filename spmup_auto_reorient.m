function spmup_auto_reorient(p)

% simple function to realign each images with the template
% only use affine transform - this is the same idea as setting the AC point
% manually except that (0,0,0) is now the same as the template (a priori
% better) - taken from the web but can't find the site again, if you know
% who wrote something like this please let me know so I can credit that
% person properly ...
% 
% Cyril Pernet June 2015
% ------------------------
% spmup copyright

if ~nargin
    [p,sts] = spm_select(Inf,'image','Select the image(s) to re-orient');
    if ~sts, return; end
end

p = cellstr(p);
vg = spm_vol(fullfile(spm('Dir'),'canonical','avg152T1.nii'));
tmp = [tempname '.nii'];

for i=1:numel(p)
    spm_smooth(p{i},tmp,[12 12 12]);
    vf = spm_vol(tmp);
    M  = spm_affreg(vg,vf,struct('regtype','rigid'));
    [u,s,v] = svd(M(1:3,1:3));
    M(1:3,1:3) = u*v';
    N  = nifti(p{i});
    N.mat = M*N.mat;
    create(N);
end
spm_unlink(tmp);
