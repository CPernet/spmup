function test_suite = test_censoring %#ok<*STOUT>
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function value = input_file()
    test_data_folder = fullfile(root_test_folder(), 'data');
    value = spm_select('FPlistRec', test_data_folder, '^rp_sub.*bold.txt');
end

function test_basic()

    data = spm_load(input_file());
    censoring_regressors = spmup_censoring(data);

    assert(size(censoring_regressors, 2) == 1);

    teardown();

end

function teardown()
    delete(fullfile(root_test_folder(), 'data', 'sub-01', 'func', '*design.txt'));
end
