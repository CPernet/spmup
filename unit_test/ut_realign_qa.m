function tests = ut_realign_qa
tests = functiontests(localfunctions);
end

function value = input_file()
test_data_folder = fullfile(test_folder(), 'data');
value = spm_select('FPlistRec', test_data_folder, '^sub.*bold.nii');
end

function value = plot_figures()
    value = 'off';
end

function test_motion(testCase) 

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
end

function test_FD(testCase) 

new_files = spmup_realign_qa(input_file(), ...
                             'Framewise displacement', 'on', ...
                             'Globals', 'off', ...
                             'Voltera', 'off', ...
                             'Movie', 'off', ...
                             'figure', plot_figures());

% 6 motion + FD + 1 censoring regressors
motion_and_fd_censor = spm_load(new_files{1});
assert(size(motion_and_fd_censor, 2) == 8);
% make sure all censoring regressors are at the end
assert(sum(motion_and_fd_censor(:, end)) == 1);

metadata = spm_load(spm_file(new_files{1}, 'ext', '.json'));
assert(numel(metadata.Columns) == 8);
end

function test_volterra(testCase) 

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
end

function test_globals(testCase) 

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
end

function test_all_together(testCase) 

new_files = spmup_realign_qa(input_file(), ...
                             'Framewise displacement', 'on', ...
                             'Globals', 'on', ...
                             'Voltera', 'on', ...
                             'Movie', 'off', ...
                             'figure', plot_figures());

% 24 voltera + RMS + FD + global + 3 censoring regressors
all_regressors = spm_load(new_files{1});
assert(size(all_regressors, 2) == 27);
% make sure all censoring regressors are at the end
assert(sum(all_regressors(:, end)) == 1);

metadata = spm_load(spm_file(new_files{1}, 'ext', '.json'));
assert(numel(metadata.Columns) == 27);

metadata = spm_load(spm_file(new_files{1}, 'ext', '.json'));
assert(strcmp(metadata.Columns{end}, 'outlier_0001'));
end


%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
this_path = fileparts(mfilename('fullpath'));
addpath(fullfile(this_path, 'utils'));
if is_github_ci()
    root_dir = getenv('GITHUB_WORKSPACE');
    addpath(fullfile(root_dir, 'spm12'));
    run(fullfile(this_path, '..'), spmup());
end
end

function setup(testCase)  % do not change function name
delete(fullfile(test_folder(), 'data', 'sub-01', 'func', '*design.txt'));
delete(fullfile(test_folder(), 'data', 'sub-01', 'func', '*.ps'));
delete(fullfile(test_folder(), 'data', 'sub-01', 'func', '*design.json'));
end

function teardown(testCase)  % do not change function name
delete(fullfile(test_folder(), 'data', 'sub-01', 'func', '*design.txt'));
delete(fullfile(test_folder(), 'data', 'sub-01', 'func', '*.ps'));
delete(fullfile(test_folder(), 'data', 'sub-01', 'func', '*design.json'));
end
