%% Main function to generate tests
function tests = exampleTest
tests = functiontests(localfunctions);
end

function test_smoke_test(testCase)

test_data_folder = fullfile(test_folder(), 'data');
time_series = spm_select('FPlistRec', test_data_folder, '^sub.*bold.nii');
plot_figures = 'off';

%%
new_files = spmup_realign_qa(time_series, ...
                             'Motion Parameters', 'off', ...
                             'Framewise displacement', 'off', ...
                             'Globals', 'off', ...
                             'Voltera', 'off', ...
                             'Movie', 'off', ...
                             'figure', plot_figures);

assert(isempty(new_files));
delete(new_files{1});
end

function test_smoke_motion(testCase) 

test_data_folder = fullfile(test_folder(), 'data');
time_series = spm_select('FPlistRec', test_data_folder, '^sub.*bold.nii');
plot_figures = 'off';

%%
new_files = spmup_realign_qa(time_series, ...
                             'Motion Parameters', 'on', ...
                             'Framewise displacement', 'off', ...
                             'Globals', 'off', ...
                             'Voltera', 'off', ...
                             'Movie', 'off', ...
                             'figure', plot_figures);

assert(isempty(new_files));
delete(new_files{1});
end

function test_smoke_FD(testCase) 

test_data_folder = fullfile(test_folder(), 'data');
time_series = spm_select('FPlistRec', test_data_folder, '^sub.*bold.nii');
plot_figures = 'off';

%%
new_files = spmup_realign_qa(time_series, ...
                             'Motion Parameters', 'off', ...
                             'Framewise displacement', 'on', ...
                             'Globals', 'off', ...
                             'Voltera', 'off', ...
                             'Movie', 'off', ...
                             'figure', plot_figures);

% 6 motion + FD + RMS + 3 censoring regressors
motion_and_fd_censor = spm_load(new_files{1});
assert(size(motion_and_fd_censor, 2) == 11);
% make sure all censoring regressors are at the end
assert(all(sum(all_regressors(:, end - 2:end)) == [1 1 1]));
delete(new_files{1});
end

%% Optional file fixtures  
function teardown(testCase)  % do not change function name
% close figure, for example
delete(fullfile(test_folder(), 'data', 'sub-01', 'func', '*.ps'));
delete(fullfile(test_folder(), 'data', 'sub-01', 'func', '*design.json'));
end