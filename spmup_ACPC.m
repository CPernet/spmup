function out = spmup_ACPC(varargin)

% simple routine to set the origin of your images
% based on auto_reorient.m from Kiyotaka Nemoto
% http://www.nemotos.net/?p=17
% 
% FORMAT spmup_ACPC
%        spmup_ACPC(P)
%        spmup_ACPC(P1,P2)
%
% INPUT empty - will ask to select images 
%       P a single image name or an array (see spm_select)
%         if P is an array, it will iterate through images
%       P1, P2 expect P1 a single image (source) and P2 an array
%         with 2 inputs, it computes how to align P1 to the template
%         and applies the transform to all images in P2 (P1 can be e.g.
%         the mean image coming out of realign)
%
% OUTPUT out is simple text telling you if the script work or not
%
% Cyril Pernet 
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2015


if isempty(nargin)
    P = spm_select(Inf,'image')
else
    
end

spmDir=which('spm');
spmDir=spmDir(1:end-5);
tmpl=[spmDir 'canonical/avg152T1.nii'];
vg=spm_vol(tmpl);
flags.regtype='rigid';
p=spm_select(inf,'image');
for i=1:size(p,1)
    f=strtrim(p(i,:));
    spm_smooth(f,'temp.nii',[12 12 12]);
    vf=spm_vol('temp.nii');
    [M,scal] = spm_affreg(vg,vf,flags);
    M3=M(1:3,1:3);
    [u s v]=svd(M3);
    M3=u*v';
    M(1:3,1:3)=M3;
    N=nifti(f);
    N.mat=M*N.mat;
    create(N);
end