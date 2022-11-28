% test bids processing pipeline on SPM's face repetition single subject dataset

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

%% get data with bids-matlab

bids_matlab_path = which('bids.layout');
if isempty(bids_matlab_path)
    if ~isfolder(fullfile(pwd, 'bids-matlab)'))
        system('git clone --branch dev --depth 1 https://github.com/bids-standard/bids-matlab.git');
    end
    bids_matlab_path = fullfile(pwd, 'bids-matlab');
    addpath(bids_matlab_path);
else
    bids_matlab_path = fileparts(fileparts(bids_matlab_path));
end

% needed to get function to bidsify the dataset
addpath(fullfile(bids_matlab_path, 'demos', 'spm', 'facerep', 'code'));

output_dir = bids.util.download_ds('source', 'spm', ...
                                   'demo', 'facerep', ...
                                   'out_path', fullfile(pwd, 'facerep', 'source'), ...
                                   'force', false, ...
                                   'verbose', true, ...
                                   'delete_previous', false);

BIDS_dir = convert_facerep_ds(fullfile(pwd, 'facerep', 'source'), ...
                              fullfile(pwd, 'facerep', 'raw'));

%%
BIDS_dir = fullfile(pwd, 'facerep', 'raw');

options = spmup_getoptions(BIDS_dir);

% set how many cores to use or don't and it uses N-1;
options.Ncores = 1;

% depends what you have, used for multispectral segmentation
options.anat = {'T1w'};
options.task = {'facerepetition'};

[BIDS, subjects] = spmup_BIDS_unpack(BIDS_dir, options);

[subjects, opt] = run_spmup_bids(BIDS, subjects, options);
