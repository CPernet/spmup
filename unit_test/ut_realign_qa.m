%% unit testing for spmup_realign_qa
% 

%% get data
test_data_folder = fullfile(test_folder(), 'data');

time_series = spm_select('FPlistRec', test_data_folder, '^sub.*bold.nii');

%% call
tSNR = spmup_realign_qa(time_series, ...
                        'Motion Parameters', 'on', ...
                        'Globals', 'off',...
                        'Voltera', 'off', ...
                        'Framewise displacement','off', ...
                        'Movie', 'off', ...
                        'figure', 'off');




