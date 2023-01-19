%
% (C) Copyright 2021 bidspm developers

root_dir = getenv('GITHUB_WORKSPACE');

addpath(fullfile(root_dir, 'spm12'));
addpath(fullfile(root_dir, 'MOcov', 'MOcov'));

cd(fullfile(root_dir, 'MOxUnit', 'MOxUnit'));
run moxunit_set_path();

cd(root_dir);
spmup;

addpath('unit_test/utils')

cd(root_dir);
logger('INFO', sprintf('Home is "%s"\n', getenv('HOME')));

spm('defaults', 'fMRI');

subfolder = '';
testFolder = fullfile(pwd, 'unit_test', subfolder);
success = moxunit_runtests(testFolder);

if success
  system('echo 0 > test_report.log');
else
  system('echo 1 > test_report.log');
end
