function tests = ut_censoring
tests = functiontests(localfunctions);
end

function value = input_file()
test_data_folder = fullfile(test_folder(), 'data');
value = spm_select('FPlistRec', test_data_folder, '^rp_sub.*bold.txt');
end

function value = plot_figures()
    value = 'off';
end

function test_basic(testCase)

data = spm_load(input_file());
[design,censoring_regressors] = spmup_censoring(data);

assert(size(design, 2) == 7)
assert(size(censoring_regressors, 2) == 1)

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

function teardown(testCase)  % do not change function name
delete(fullfile(test_folder(), 'data', 'sub-01', 'func', '*design.txt'));
end