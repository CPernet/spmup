function first_level_qa

% routine calling realign QA and normalize QA
% peforms the 'usual' stuff ie 
% --> plot motion paramters and 1st derivatives
% --> compute and plot globals (see spm_globals)
% --> generate new regressor.txt file with motion parameters and outliers detected using Carling's modified boxplot rule
% --> create the average of normalized images
% --> compute the difference between the T1 (assumed normalized) and the normalized images - plot the outlines
% --> make a movie in all directions assuming [0 0 0] using the template

spm('Defaults','fmri')
current = pwd;

% images as input
Normalized = spm_select(Inf,'image','select normalized images');
[folder,name,ext] = fileparts(Normalized(1,:));
cd(folder); cd ..
T1 = spm_select(1,'image','select normalized T1 image');
cd(current)

% flags
flags.motion_parameters = 'on'; 
flags.globals = 'on';
flags.volume_distance = 'off';
flags.movie = 'off';
flags.AC = [46 64 37]; 
flags.average = 'on'; 
flags.T1 = 'on';

% motion test
realign_qa(Normalized,flags)

% normalization test
normalize_qa(flags,Normalized,T1)

