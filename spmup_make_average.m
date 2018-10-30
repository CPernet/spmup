function spmup_make_average

% simple routine to create an average
% Cyril Pernet 10-07-2013

Defaults = spm_get_defaults;
[t,sts] = spm_select(Inf,'image','select images');
if sts == 0
    return
else
    V = spm_vol(t);
    check = spm_check_orientations(V);
end

images = spm_read_vols(V);
average = nanmean(images,4);
newV = spm_create_vol(V(1));
[path,name,ext]=fileparts(newV.fname);
newV.fname = [path filesep 'average' ext];
newV.descrip = 'Average file create bt spm_make_average';
V = spm_write_vol(newV,average);

