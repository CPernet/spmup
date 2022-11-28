%% get data
test_data_folder = fullfile(test_folder(), 'ds000117');

subjects = fullfile(test_data_folder, 'spmup_subjects_task-facerecognition.mat');

table_name = spmup_BIDS_QCtables(subjects, 'anat');