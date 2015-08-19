function adaptive_thresholding(stat_filename, mask_filename, spm_mat_file, con_index)

% read the data (spmT and mask), call the Gamma-Gaussian Mixture code and
% then apply this threshold to save the threasholded image with topo FDR at
% q=0.05
%
% FORMAT: adaptive_thresholding(stat_filename, mask_filename, spm_mat_file, con_index)
%
% INPUT: stat_filename: full name (ie with path) of the spmT image
%        mask_filename: full name (ie with path) of the mask image
%        spm_mat_file: full name (ie with path) of the SPM.mat
%        con_index: the number of the contrast image
% 
% see also ggmm_thresholding topo_fdr_thresholding
% ------------------------------------------
% Chris Gorgolewski 29 June 2012

% check inputs
  if nargin == 0
      stat_filename = spm_select(1, 'image', 'Select spmT-map', [], pwd, 'spmT*.*');
  end;
  
  if nargin <= 1
      % deriving other information from first input
      [p, nm, e, ~] = spm_fileparts(stat_filename);
      mask_filename = [p filesep 'mask' e];
      spm_mat_file  = [p filesep 'SPM.mat'];
      con_index     = str2num(nm(end-3:end));
  end
  
[thr,bic] = ggmm_thresholding(stat_filename, mask_filename);
if length(thr) == 1
    topo_fdr_thresholding(spm_mat_file, con_index, thr, 'none', 1, 1, 0.05, stat_filename, 'stat', 0);
elseif length(thr) == 2
    topo_fdr_thresholding(spm_mat_file, con_index, thr(2), 'none', 1, 1, 0.05, stat_filename, 'stat', 0);
elseif isempty(thr)
    disp('no threshold found ! somethings strange in the data, please check')
end
