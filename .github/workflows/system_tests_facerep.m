%
% (C) Copyright 2021 bidspm developers

root_dir = getenv('GITHUB_WORKSPACE');

fprintf('\nroot dir is %s\n', root_dir);

addpath(fullfile(root_dir, 'spm12'));

cd(fullfile(root_dir, 'system_test'));

run facerep;
