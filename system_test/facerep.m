% test bids processing pipeline on SPM's face repetition single subject dataset

%% update path
local_test_dir = '/indirect/data1/cpernet/';

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
if ~exist(fullfile(bids_matlab_path, 'demos', 'spm', 'facerep', 'code'),'dir')
    command = ['cd ' bids_matlab_path ' && git switch dev'];
    system(command)
end        
addpath(fullfile(bids_matlab_path, 'demos', 'spm', 'facerep', 'code'));

try
    if isempty(local_test_dir)
        output_dir = bids.util.download_ds('source', 'spm', ...
            'demo', 'facerep', ...
            'out_path', fullfile(pwd, 'facerep', 'source'), ...
            'force', false, ...
            'verbose', true, ...
            'delete_previous', true);
    else
        output_dir = bids.util.download_ds('source', 'spm', ...
            'demo', 'facerep', ...
            'out_path', local_test_dir, ...
            'force', false, ...
            'verbose', true, ...
            'delete_previous', true);
    end
    BIDS_dir = convert_facerep_ds(fullfile(pwd, 'facerep', 'source'), ...
        fullfile(pwd, 'facerep', 'raw'));
    
catch downloaderr
    if ~isempty(local_test_dir)
        if exist(fullfile(local_test_dir,'face_rep'),'dir')
            output_dir = fullfile(local_test_dir,'face_rep');
            warning('\ncould not download the facerep data automatically, %s %s',downloaderr.identifier,downloaderr.message)
            warning('using found directory %s',output_dir)
            BIDS_dir = convert_facerep_ds(output_dir, fullfile(output_dir, 'raw'));
        else
            error('\ncould not download the facerep data automatically, %s %s',downloaderr.identifier,downloaderr.message)
        end
    end
end


%% now run spmup on the face_rep data

options            = spmup_getoptions(BIDS_dir);
options.Ncores     = 1;% set how many cores to use or don't and it uses N-1;
options.anat       = {'T1w'};
options.task       = {'facerepetition'};
options.conditions = 'face_type'; % that's the conditions to model
[BIDS, subjects]   = spmup_BIDS_unpack(BIDS_dir, options);
[subjects, opt]    = run_spmup_bids(BIDS, subjects, options);
