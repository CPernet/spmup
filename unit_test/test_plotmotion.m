function test_suite = test_plotmotion %#ok<*STOUT>
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function test_smoke()

    test_data_folder = fullfile(root_test_folder(), 'data');
    input_file.name = spm_select('FPlistRec', test_data_folder, '^rp_sub.*bold.txt');

    spmup_plotmotion(input_file);

end
