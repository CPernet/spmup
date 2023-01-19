function test_suite = test_realign_qa %#ok<*STOUT>
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function value = input_file()
    test_data_folder = fullfile(root_test_folder(), 'data');
    value = spm_select('FPlistRec', test_data_folder, '^sub.*bold.nii');
end

function value = plot_figures()
    value = 'off';
end

function test_motion()

    teardown();

    new_files = spmup_realign_qa(input_file(), ...
                                 'Framewise displacement', 'off', ...
                                 'Globals', 'off', ...
                                 'Voltera', 'off', ...
                                 'Movie', 'off', ...
                                 'figure', plot_figures());

    motion_parameters = spm_load(new_files{1});
    assert(size(motion_parameters, 2) == 6);

    metadata = spm_load(spm_file(new_files{1}, 'ext', '.json'));
    assert(numel(metadata.Columns) == 6);

    teardown();

end

function test_FD()

    teardown();

    new_files = spmup_realign_qa(input_file(), ...
                                 'Framewise displacement', 'on', ...
                                 'Globals', 'off', ...
                                 'Voltera', 'off', ...
                                 'Movie', 'off', ...
                                 'figure', plot_figures());

    % 6 motion + FD + RMS + 3 censoring regressors
    motion_and_fd_censor = spm_load(new_files{1});
    assert(size(motion_and_fd_censor, 2) == 11);
    % make sure all censoring regressors are at the end
    assert(all(sum(motion_and_fd_censor(:, end - 2:end)) == [1 1 1]));

    metadata = spm_load(spm_file(new_files{1}, 'ext', '.json'));
    assert(numel(metadata.Columns) == 11);

    teardown();

end

function test_volterra()

    teardown();

    new_files = spmup_realign_qa(input_file(), ...
                                 'Framewise displacement', 'off', ...
                                 'Globals', 'off', ...
                                 'Voltera', 'on', ...
                                 'Movie', 'off', ...
                                 'figure', plot_figures());

    % 6 motion + their derivatives + square of each
    voltera = spm_load(new_files{1});
    assert(size(voltera, 2) == 24);

    metadata = spm_load(spm_file(new_files{1}, 'ext', '.json'));
    assert(numel(metadata.Columns) == 24);

    teardown();

end

function test_globals()

    teardown();

    new_files = spmup_realign_qa(input_file(), ...
                                 'Framewise displacement', 'off', ...
                                 'Globals', 'on', ...
                                 'Voltera', 'off', ...
                                 'Movie', 'off', ...
                                 'figure', plot_figures());

    % 6 motion + one global regressor
    globals = spm_load(new_files{1});
    assert(size(globals, 2) == 7);

    metadata = spm_load(spm_file(new_files{1}, 'ext', '.json'));
    assert(numel(metadata.Columns) == 7);

    teardown();

end

function test_all_together()

    teardown();

    new_files = spmup_realign_qa(input_file(), ...
                                 'Framewise displacement', 'on', ...
                                 'Globals', 'on', ...
                                 'Voltera', 'on', ...
                                 'Movie', 'off', ...
                                 'figure', plot_figures());

    % 24 voltera + + FD  + RMS  + global + 3 censoring regressors
    all_regressors = spm_load(new_files{1});
    assert(size(all_regressors, 2) == 30);
    % make sure all censoring regressors are at the end
    assert(all(sum(all_regressors(:, end - 2:end)) == [1 1 1]));

    metadata = spm_load(spm_file(new_files{1}, 'ext', '.json'));
    assert(numel(metadata.Columns) == 30);

    metadata = spm_load(spm_file(new_files{1}, 'ext', '.json'));
    assert(strcmp(metadata.Columns{end}, 'outlier_0003'));

    teardown();

end

function teardown()
    delete(fullfile(root_test_folder(), 'data', 'sub-01', 'func', '*design.txt'));
    delete(fullfile(root_test_folder(), 'data', 'sub-01', 'func', '*.ps'));
    delete(fullfile(root_test_folder(), 'data', 'sub-01', 'func', '*design.json'));
end
