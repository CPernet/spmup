% test bids processing pipeline on the ds000117 dataset

%% update path

% in case spmup was not added to the path
location = fullfile(pwd, '..');
addpath(fullfile(location, 'adaptative_threshold'), ...
        fullfile(location, 'bids'), ...
        fullfile(location, 'external'), ...
        fullfile(location, 'hrf'), ...
        fullfile(location, 'plot'), ...
        fullfile(location, 'QA'), ...
        fullfile(location, 'utilities'));

BIDS_dir = fullfile(pwd, 'ds000117');

%%
options = spmup_getoptions(BIDS_dir);

% set how many cores to use or don't and it uses N-1;
options.Ncores = 6;

% depends what you have, used for multispectral segmentation
options.anat = {'T1w'};
options.task = {'facerecognition'};
options.subjects = {'sub-01', 'sub-02'};

[BIDS, subjects] = spmup_BIDS_unpack(BIDS_dir, options);

[subjects, opt] = run_spmup_bids(BIDS, subjects, options);
