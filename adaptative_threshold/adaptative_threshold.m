function T = adaptative_threshold(varargin)

% evaluate from an SPM T map the mixing distributions corresponding to the 
% null hypothesis (in theory a Gaussian centred on 0, but the algorithms 
% allows a shift because of global effects) and the hypotheses of activation
% and deactivation (a Gamma and a negated Gamma). Once the data fitting is
% done, the uncorrected p value to be used in spm/results is returned -
% results can be then evaluated based on the FDR extend cluster threshold
% 
% FORMAT: T = adaptative_threhold(P,SPM)
%
% INPUT: - this can be left empty and the user is prompted to pick an image,
%        and the flag will be set to 1
%        - P is the full name (ie with the path) of the SPM T image
%        - SPM is the SPM.mat corresponding to the image P, if not
%        specified one looks for an SPM in the same folder as the T image
%
% OUTPUT: - T is the t threshold corresponding to the crossing point between
%         the Gaussian and the positive Gamma
%
% This code is an adaptation of the python code available with Nipype.
% Reference: Gorgolewski et al. Adaptive thresholding for reliable 
% topological inference in single subject fMRI analysis - Front. Hum. 
% Neurosci. 6:245. doi: 10.3389/fnhum.2012.00245
%
% see also ggmm_thresholding topo_fdr_thresholding
% ------------------------------------------------
% Cyril Pernet and Chris Gorgolewski v1 - July 2012


% This file is part of the adaptative_threshold SPM plug-in
%
% adaptative_threshold SPM plug-in is free software: you can redistribute 
% it and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation, either version 3 of the 
% License, or (at your option) any later version.
%
% adaptative_threshold SPM plug-in is distributed in the hope that it will 
% be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with adaptative_threshold SPM plug-in.  
%     If not, see <http://www.gnu.org/licenses/>.
%
% -----------------------------
% Copyright (C) Cyril Pernet and Chris Gorgolewski 2012

disp(' ')
disp('This file is part of the adaptative_threshold SPM plug-in')
disp(' ')
disp('adaptative_threshold SPM plug-in is free software: you can redistribute it')
disp('and/or modify it under the terms of the GNU General Public License as')
disp('published by the Free Software Foundation, either version 3 of the')
disp('License, or (at your option) any later version.')
disp(' ')
disp('adaptative_threshold SPM plug-in is distributed in the hope that it will ')
disp('be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of')
disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.')
disp(' ')
disp('See the GNU General Public License for more details.')
disp('You should have received a copy of the GNU General Public License')
disp('along with adaptative_threshold SPM plug-in.  ')
disp('If not, see <http://www.gnu.org/licenses/>.')
disp(' ')



%% check inputs/outputs
% ---------------------
SPM =[];
if isempty(varargin)
    [P,sts] = spm_select(1,'image','Select spmT image');
    if sts == 0
        return
    end
elseif nargin == 1
    P = varargin{1};
elseif nargin == 2
    P = varargin{1};
    SPM = varargin{2};
elseif nargin > 2
    error('too many input arguments')
end

% get the SPM.mat
if isempty(SPM)
    try
        SPM_path = P(1:max(find(P=='/')));
        cd(SPM_path)
        load('SPM.mat')
    catch ME
        disp('SPM.mat cannot be located')
        [S,sts] = spm_select(1,'^SPM.*\.mat$','Select SPM.mat');
        if sts == 0
            return
        else
            load(S);
            SPM_path = SPM.swd;
        end
    end
end


% check for mask
mask_choice = spm_input('Masking..','!+1','b',{'load mask','implicit'});
if strcmp(mask_choice,'implicit')
    mask_name = [pwd '/mask.nii'];
    try
        mask_data = spm_read_vols(spm_vol(mask_name));
    catch me
        [mask_name,sts] = spm_select(1,'image','Coudn''t locate the default mask, please select mask image');
        if sts == 0
            return
        else
            mask_data = spm_read_vols(spm_vol(mask_name));
        end
    end
else
    [mask_name,sts] = spm_select(1,'image','Select mask image');
    if sts == 0
        return
    else
        mask_data = spm_read_vols(spm_vol(mask_name));
    end
end

% get some info on this image
% ----------------------------
V = spm_vol(P);
data = spm_read_vols(V);
data = data(mask_data > 0);

% get the information about this contrast
[pth,nam,ext,num] = spm_fileparts(P);
contrast_nb = eval(nam(6:end));
if ~strcmp(SPM.xCon(contrast_nb).STAT,'T')
    error('mismatch between SPM.mat and spmT image'); % not a T stat
end


%% fit the Gamma-Gaussian-Gamma mixture model and make a figure
% ------------------------------------------

[threshold,bic] = ggmm_thresholding(P, mask_name);
[value,winner] = min(bic);


%% return the outputs
% ---------------------

if (winner == 2 || winner == 3) && ~isempty(threshold)
    % df = SPM.xCon(contrast_nb).eidf; % always 1 for T
    dfe = SPM.xX.erdf;
    % p = 1-tcdf(threshold,dfe);
    if nargout == 0
        msgbox(['For cluster inference using FDR you can set the T threashold @ ' num2str(threshold) ' (set p value adjustment to no) - only clusters with corrected p values must be used'],'Cluter forming threshold','help')
    else
        T = threshold;
    end
else
    msgbox('there is no evidence for a positive signal','Cluter forming threshold','warn')
    T = NaN;
end


%% additional feature - save a threasholded image
% ----------------------------------------
    
saving_choice = spm_input('Saving the image threasholded with topological FDR?','!+1','b',{'yes','no'});
if strcmp(saving_choice,'yes')
    topo_fdr_thresholding([SPM_path 'SPM.mat'], contrast_nb, threshold(end), 'none', 1, 1, 0.05, P, 'stat', 0);
end


