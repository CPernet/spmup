%% unit testing for spmup_temporalSNR
% 
% the model is simple with 4 'boxes': GM, WM, CSF, background
% noise is white noise everywhere (std = 1), and the signal is a simple cos
% added with amplitude 1 (background) 2 (GM) 3 (WM) and 4 (CSF)
% ------------------------------
% Cyril Pernet 05 January 2016
   
%% make up the data

% pick up SPM template and populate



%% call
tSNR = spmup_temporalSNR(time_series,masks,0);

%% check we have what we expect








