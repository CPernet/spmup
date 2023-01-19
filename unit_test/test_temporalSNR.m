function test_suite = test_temporalSNR %#ok<*STOUT>
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function test_suite = test_smoke

    %% unit testing for spmup_temporalSNR
    %
    % the model is simple with 4 'boxes': GM, WM, CSF, background
    % noise is white noise everywhere (std = 1), and the signal is a simple cos
    % added with amplitude 1 (background) 2 (GM) 3 (WM) and 4 (CSF)
    % ------------------------------
    % Cyril Pernet 05 January 2016

    %% get data
    test_data_folder = fullfile(root_test_folder(), 'data');

    time_series = spm_select('FPlistRec', test_data_folder, '^sub.*bold.nii');
    gm = spm_select('FPlistRec', test_data_folder, '^rc1.*.nii');
    wm = spm_select('FPlistRec', test_data_folder, '^rc2.*.nii');
    csf = spm_select('FPlistRec', test_data_folder, '^rc3.*.nii');

    %% call
    tSNR = spmup_temporalSNR(time_series, ...
                             {gm; wm; csf}, ...
                             'figure', 'save', ...
                             'SNR0', 'on', ...
                             'linearity', 'on');

end
