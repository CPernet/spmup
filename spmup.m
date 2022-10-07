function spmup

% toolbox setting - run it once to set all paths
% ---------------------------------------------
% Copyright (c) SPM Utility Plus toolbox

location = fileparts(which('spmup.m'));
addpath(fullfile(location,'adaptative_threshold'), ...
    fullfile(location,'bids'), ...
    fullfile(location,'external'), ...
    fullfile(location,'hrf'), ...
    fullfile(location,'plot'), ...
    fullfile(location,'QA'), ...
    fullfile(location,'utilities'));

if ~exist(fullfile(fileparts(location),'gp_event_plot'),'dir')
    mkdir(fullfile(fileparts(location),'gp_event_plot'))
    copyfile(fullfile([location filesep 'plot'],'gp_event_plot.m'),...
        fullfile(fileparts(location),'gp_event_plot'))
end
savepath

disp('----------------------------------------------------')
disp('            SPM Ulitity Plus toolbox                ')
disp('----------------------------------------------------')
disp(' standard options are available using the SPM batch ')
disp(' for the pipeline see run_spmup_bids ')
help run_spmup_bids
disp('----------------------------------------------------')
