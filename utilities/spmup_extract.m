function timeseries = spmup_extract(fmridata,gaymatter,atlas)

% routine to extract timeseries data as per atlas
%
% FORMAT timeseries = spmup_extract(fmridata,gaymatter,atlas)
%
% INPUT
%
% OUTPUT timeseries is a cell array of voxels for each ROI of the atlas intersecting
%                   with gray matter, dimensionis X voxels * Y volumes
%
% Cyril Pernet 
% --------------------------
%  Copyright (C) SPMUP Team 

if nargin == 0
    
    [fmridata, sts] = spm_select(1, 'image', 'Select fMRI data');
    if sts == 0
        disp('selection aborded');
        return
    end
    if size(fmridata, 1) == 1 && strcmp(fmridata(length(fmridata) - 1:end), ',1')
        fmridata = fmridata(1:length(fmridata) - 2); % in case picked 4D put left ,1
    end
    
    [cn1, sts] = spm_select(1, 'image', 'Select grey matter tissue image');
    if sts == 0
        disp('selection aborted');
        return
    end
    
    [cn2, sts] = spm_select(1, 'image', 'Select atlas image image');
    if sts == 0
        disp('selection aborted');
        return
    end
end

% check fmridata is 4D
if iscell(fmridata)
    
    if strcmp(fmridata{1}(end - 1:end), ',1')
        fmridata{1} = fmridata{1}(1:end - 2);
    end
    
    if strcmp(fmridata{2}(end - 1:end), ',1')
        fmridata{2} = fmridata{2}(1:end - 2);
    end
    
    % get the data
    Vfmri = spm_vol(fmridata{1});
    
else % not a cell array
    if strcmp(fmridata(1, end - 1:end), ',1')
        fmridata = fmridata(1:end - 2);
    end
    Vfmri = spm_vol(fmridata(1, :));
end

filepath = fileparts(Vfmri(1).fname);
disp   ('------------------------------')
fprintf(' running spmup_extract on %s\n',Vfmri(1).fname)
disp   ('------------------------------')

%% deal with the masks
if iscell(cn1)
    cn1 = cell2mat(cn1);
end

greymatter = spm_vol(cn1);
if any(greymatter.dim ~= Vfmri(1).dim)
    error('Dimension issue between data and grey matter');
end

c1 = spm_read_vols(spm_vol(greymatter));
if min(c1(:)) < 0.5
    warning('the grey matter mask has values as low as %g, you may consider thresholding it',min(c1(:)))
end

if iscell(cn2)
    cn2 = cell2mat(cn2);
end

if any(cn2.dim ~= Vfmri(1).dim)
    warning('Dimension issue between data and atlas - interpolating the atlas');
    
    clear matlabbatch;
    matlabbatch{1}.spm.util.bbox.image = {Vfmri(1).fname};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    bb = spm_jobman('run', matlabbatch);
    copyfile(cn2.fname,fullfile(fileparts(Vfmri(1).fname), atlas_filename));
    roi = fullfile(fileparts(Vfmri(1).fname), atlas_filename);
    Vroi = spm_vol(cell2mat( ...
        spmup_resize(roi, bb{1}.bb, ...
        abs(diag(Vfmri(1).mat(1:3, 1:3)))')));
end
atlas = spm_read_vols(spm_vol(c2));


