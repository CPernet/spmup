function spmup_find_dropout(varargin)

% Simple routine to loop throughout masks - create the sum of masks rather 
% than the intersection, search through slices (z direction) to see if some 
% subjects are outliers, and look at total overlap to see again if some 
% subjects are outliers - the focus is on Z direction as it is often easier
% to spot normalization issue this way
%
% OUTPUTS: sum_of_masks.img 
%          --> the sum of masks to be compared with the actual mask which
%          is the intersection of masks
%          --> mask_drop_outz.eps the figure diplays which slices are
%          included in the final mask - allowing to spot bad subjects
%          --> in the command window, also report bad subjects identified
%          using a box plot rule as subjects whith many more slices not
%          included in the final mask
%
% Inpited by find_dropout by Jason Steffener
%
% Cyril Pernet 13-03-13

%% get the data
% get the 2nd level analysis
if nargin == 0
[fname,pname]=uigetfile('*SPM.mat','Select 2nd level design matrix');
else
[pname,fname,ext]=fileparts(varargin);
fname = [fname ext];
end
cd(pname); current = pwd;
V = spm_vol('mask.nii');
mask = spm_read_vols(V);
load(fullfile(pname,fname));

% go through each volumes to make sure we have all subjects
for subject = 1:SPM.nscan
    [folder{subject},name,ext] = fileparts(SPM.xY.P{subject});
end
folder = unique(folder); % retain unique folders

% loop through folders to load the masks of the 1st level
N=size(folder,2);
data = zeros([V.dim N]);
for subject=1:N
    cd(folder{subject});
    try
        data(:,:,:,subject) = spm_read_vols(spm_vol('mask.nii'));
    catch
        cd ..
        data(:,:,:,subject) = spm_read_vols(spm_vol('mask.nii'));
    end
end

%% sum masks 
cd(current); mkdir('find_dropout_output'); cd('find_dropout_output')
final = sum(data,4);
V.fname = [pwd filesep 'sum_of_masks.img'];
spm_write_vol(V,final);

%% slice by slice find differences (in z direction)
record = zeros(V.dim(3),N);
for v = 1:V.dim(3)
    current_data = squeeze(final(:,:,v));
    number_subjects(v) = max(current_data(:));
    if number_subjects(v) ~=N % which subject is not participating
        for subject = 1:N
            if isempty(find(squeeze(data(:,:,v,subject))))
                record(v,subject) = 1;
            end
        end
    end
end

% check compared to the group if some subject stand out
% record = nb of slices present per subject
bad_subjects = detect_outliers(sum(record)'); % tells if a subject has doesn't contribute to the mask
if isempty(bad_subjects)
    disp('all subjects show similar slice overalp in z direction')
else
    for s=1:length(bad_subjects)
        fprintf('compared to the group, subject %g \n %s \n is noticeably badly normalized in z \n',bad_subjects(s), folder{bad_subjects(s)});
    end
end
disp(' ')


figure; subplot(1,5,1:4);
imagesc(record); colormap('gray')
xlabel('subjects','Fontsize',14)    
ylabel('slices','Fontsize',14)    
title('bad slices (z) per subject','Fontsize',14)    
subplot(1,5,5); tmp = (sum(squeeze(sum(mask,1)),1))';
imagesc(tmp == 0)   
title('Current mask','Fontsize',14)    
saveas(gcf, 'mask_drop_outz.eps','psc2'); 


%% slice by slice find differences (in y direction)
record = zeros(V.dim(2),N);
for v = 1:V.dim(2)
    current_data = squeeze(final(:,v,:));
    number_subjects(v) = max(current_data(:));
    if number_subjects(v) ~=N % which subject is not participating
        for subject = 1:N
            if isempty(find(squeeze(data(:,v,:,subject))))
                record(v,subject) = 1;
            end
        end
    end
end

% check compared to the group if some subject stand out
% record = nb of slices present per subject
bad_subjects = detect_outliers(sum(record)'); % tells if a subject has doesn't contribute to the mask
if isempty(bad_subjects)
    disp('all subjects show similar slice overalp in y direction')
else
    for s=1:length(bad_subjects)
        fprintf('compared to the group, subject %g \n %s \n is noticeably badly normalized in y \n',bad_subjects(s), folder{bad_subjects(s)});
    end
end
disp(' ')

figure; subplot(1,5,1:4);
imagesc(record); colormap('gray')
xlabel('subjects','Fontsize',14)    
ylabel('slices','Fontsize',14)    
title('bad slices (y) per subject','Fontsize',14)    
subplot(1,5,5); tmp = sum(squeeze(sum(mask,1)),2);
imagesc(tmp == 0)   
title('Current mask','Fontsize',14)    
saveas(gcf, 'mask_drop_outy.eps','psc2'); 

%% slice by slice find differences (in x direction)
record = zeros(V.dim(1),N);
for v = 1:V.dim(1)
    current_data = squeeze(final(v,:,:));
    number_subjects(v) = max(current_data(:));
    if number_subjects(v) ~=N % which subject is not participating
        for subject = 1:N
            if isempty(find(squeeze(data(v,:,:,subject))))
                record(v,subject) = 1;
            end
        end
    end
end

% check compared to the group if some subject stand out
% record = nb of slices present per subject
bad_subjects = detect_outliers(sum(record)'); % tells if a subject has doesn't contribute to the mask
if isempty(bad_subjects)
    disp('all subjects show similar slice overalp in x direction')
else
    for s=1:length(bad_subjects)
        fprintf('compared to the group, subject %g \n %s \n is noticeably badly normalized in x \n',bad_subjects(s), folder{bad_subjects(s)});
    end
end
disp(' ')

figure; subplot(1,5,1:4);
imagesc(record); colormap('gray')
xlabel('subjects','Fontsize',14)    
ylabel('slices','Fontsize',14)    
title('bad slices (x) per subject','Fontsize',14)    
subplot(1,5,5); tmp = sum(squeeze(sum(mask,2)),2);
imagesc(tmp == 0)    
title('Current mask','Fontsize',14)    
saveas(gcf, 'mask_drop_outx.eps','psc2'); 


%% volume by volume find differences in total overalp

non_zeros = numel(find(final>0)); % how many voxels with at leat one subject
for subject = 1:N
    A = squeeze(data(:,:,:,subject)) == mask;
    C(subject) =  sum(A(:)) / numel(mask) * 100;
end

bad_subjects = detect_outliers(C);
for s=1:length(bad_subjects)
        fprintf('compared to the group, subject %g \n %s \n has a very good overlap, maybe biasing the mask \n',bad_subjects(s), folder{bad_subjects(s)});
end
disp(' '); disp('Percentages of overlap = '); disp(num2str(C')); disp(' ')
disp('find dropout done')

end

function bad_subjects = detect_outliers(out)
% use the MAD which detects even small deviations

N = length(out);
M = nanmedian(out);
MAD= nanmedian(abs(out-M));
if N == 2
    bn=1.197;
elseif N == 3
    bn=1.49;
elseif N == 4
    bn=1.36;
elseif N == 5
    bn=1.217;
elseif N == 6
    bn=1.189;
elseif N == 7
    bn=1.138;
elseif N == 8
    bn=1.127;
elseif N == 9
    bn=1.101;
else
    bn=N/(N-0.8);
end
k = 2.2414; % = sqrt(chi2inv(0.975,1))      
MADS=MAD.*1.4826.*bn;
bad_subjects = find(out > M+k.*MADS);

end    
    
    
    
    
