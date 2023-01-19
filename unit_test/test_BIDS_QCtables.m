function test_suite = test_BIDS_QCtables %#ok<*STOUT>
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function test_smoke()

    test_data_folder = fullfile(root_test_folder(), 'ds000117');

    subjects = fullfile(test_data_folder, 'spmup_subjects_task-facerecognition.mat');

    table_name = spmup_BIDS_QCtables(subjects, 'anat');

end
