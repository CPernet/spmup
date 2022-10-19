# SPM UP

SPM utlities plus are tools designed to get the best of mass univariate analyses.

## Installation

This acts as a toolbox - place it in /spm12/toolbox/ and run spmup.m which will add all the other paths, copy gp_event_plot.m to appear in the spm toolbox tab, and save this.

## Pipeline

There is an automated pipeline that makes the most of SPMUP. Assuming your data are in BIDS you can use the functions located in the `bids` folder. A typical usage would follow these steps:

```matlab
BIDS_dir        = 'mybidsdatasetpath';
options         = spmup_getoptions(BIDS_dir);
options.Ncores  = N; % set how many cores to use or don't and it uses N-1;
options.anat    = {'T1w','T2w'}; % depends what you have, used for multisprectral segmentation 
[BIDS,subjects] = spmup_BIDS_unpack(BIDS_dir,options);
system(['chmod -Rf 755 ' options.outdir]) % on servers you often need to do that as matlab screw up permissions
[subjects,opt]  = run_spmup_bids(BIDS,subjects,options); % that's the magic bit
```

## QA

The QA folder contains a series of tools to check the quality of your images. Some metrics can be used at the group level as covariates, which can be important when comparing different group of subjects, if their quality metrics are different. The intented usage is to run fMRI data preprocessing and then call spmup_anatQA.m on the coregistered T1 image (i.e. the one used to derive normalization parameters), spmup_temporalSNR.m and spmup_fisrt_level_qa.m on the smoothed normalized realigned the slice timed fMRI data (i.e. the data used for your GLM). Post GLM estimation, spmup_glmQA is called. At the group level, call spmup_second_level_qa.m.

### subject level QA: SNR and beyond

_spmup_anatQA:_

Inspired by the [Preprocessed Connectome Project Quality Assurance Protocol](http://preprocessed-connectomes-project.org/quality-assessment-protocol/), this function takes an anatomical image along with gray (c1) and white (c2) matter images to returns
- SNR (Signal-to-Noise Ratio, the mean intensity within gray and white matter divided by the standard deviation of the values outside the brain),
- CNR (Contrast to Noise Ratio, the mean of the white matter intensity values minus the mean of the gray matter intensity values divided by the standard deviation of the values outside the brain),
- FBER (Foreground to Background Energy Ratio, the variance of voxels in grey and white matter divided by the variance of voxels outside the brain),
- EFC (Entropy Focus Criterion, the entropy of voxel intensities proportional to the maximum possibly entropy for a similarly sized image, indicates ghosting and head motion-induced blurring),
- an attempt on a rough asymetry measure (Left-Right/Left+Right Assuming the [0,0,0] coordinate is roughly in the middle of the image - large difference can indicate an issue with the image (for healthy brains).

_spmup_temporalSNR:_

This function recapitulates tSNR as described in [Thomas Liu (2016)](https://www.sciencedirect.com/science/article/abs/pii/S1053811916304694) and the physio to thermal noise ratio and correlation discussed in [Wald and Polimeni (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5483395/). The function is retuning:
- tSNR.GM: mean GM signal / std over time (estimate BOLD from GM>(WM+CSF))
- tSNR.WM:  mean WM signal / std over time (estimate non-BOLD from WM>(GM+CSF))
- tSNR.CSF: mean CSF signal / std over time (estimate non-BOLD from CSF>(GM+WM))
- tSNR.Background:  mean signal outside mask (GM+WM+CSF) / std over time + report the data as this should only be termal noise, i.e. gaussian distributed (figure print)
- tSNR.average (tSNR): mean signal / sqrt(std(GM)^2+std(WM+CSF)^2+std(Background)^2)
- tSNR.image (SNR0):  mean signal inside mask / std outside mask over time
- tSNR.physio2termal_ratio: sqrt((tSNR(whole image)/SNR0(brain only))^2-1)
- tSNR.physio2termal_corr: correlation between images
- tSNR.roi: tSNR for increased ROI (from in mask by increasing slices) ~linear function of srqrt(nb voxels)
- tSNR.signal_mean: sqrt(std(GM)^2+std(WM+CSF)^2) / sqrt((SNR0^2/tSNR- 1)/SNR0^2)
- a tSMR_time_series.nii image is also saved on the drive, showing tSNR in each voxel for GM, WM and CSF as computed above
- a post-scrip figure with plots of noise, tSNR and ratio with SNR0, noise per increasing ROI size (should be a line)

_spmup_spectral_:

This function computes the power spectrum slice by slice and looks for outliers across volumes. This is not intended to find time series outliers (see below) but to check for artefacts in the data (usually a coil a reconstruction issue giving rise to large spectral changes).

_spmup_sfs_:

This function computes the Signal Fluctuation Sensitivity metric [DeDora et 2016](http://journal.frontiersin.org/article/10.3389/fnins.2016.00180/full) which is arguably or more sensitive metric than tSNR to mesure BOLD variations in resting state fMRI. By default sfs is returned in the 7 resting state networks from [Yeo et al 2011](https://www.physiology.org/doi/full/10.1152/jn.00338.2011).

### fMRI 1st level QA

_spmup_fisrt_level_qa_

This high level functions returns information on files generated by other functions:
- spm_basics: mean and std images
- spmup_realign_QA: makes plots of motions and globals, creates an augmented design matrix (calling spmup_censoring) and makes movies

_spmup_censoring:_

This function takes motion parameters along with any valued vectors to create a design.txt file to be used as regressors. Along with motion parameters (and Voltrera expansion if specified), a series of binary vectors are added for outlying data. Values vectors would typically come from spmup_FD and spm_globals. Note that outliers will differ from the 'Time series outliers' functions detailled below, because spmup_censoring outliers are defined with high specificity but low sensitivity using a box-plot rule (interquartile range) with Carling's k adjustment i.e. less outliers are found but they are for sure outliers. The default of those functions uses S-outliers witch is a method with lower specificity but higher sensitivity (= more false positives, but good for QA the data).

### fMRI 2nd level QA

_spmup_second_level_qa_

This function takes a group level SPM.mat to then check for signal droupout from individual masks (spmup_find_dropout)

### Time series outliers

There are many ways to find outlying volumes in fMRI, either based on the volumes themselves or based on derived motion parameters. These are particularly useful to measure/see the effect of preprocessing. Each of these function returms vector values and logical vectors indicating outliers. Vector values can be used in spmup_censoring with uses a less sensitive metric to define outiers. Here, outliers are defined with high sensitivity but low specificity using Carling's K i.e. some declared outliers might not actually have artefacts (see spmup_comp_robust_outliers) but this is whart is needed to diagnose problems. These functions are intented to measure/see the effect of preprocessing (i.e. run on raw then on reprocessed).

_spmup_voxel_outliers:_

This function is similar to [3dToutcount](http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dToutcount.html) retuning which volumes are outliers based on the number of outlying voxels (defined as time points deviating from the median).

_spmup_volumecorr:_

This function computes the volume-wise correlations and also returns outlying volumes. Data are expected to be correlated in time, a large shift in the overall correlation indicates an issue.

_spmup_spatialcorr:_

This function compute the slice-wise correlations. This is performed both between volumes (thus similar to spmup_volumecorr) and within volumes. This returns a figure useful to spot spatial dropouts.

_spmup_FD:_

This function computes the framewise displacement from motion parameters defined as the L1 and L2 norms from derivatives, and returns outliers.

## hrf (re)estimation

Using the GLM with a temporal derivative, we can compute a systematic time shift per condition. If indeed at each trial, the reponse peaks earlier or later, then the hrf regressor does not fit the maximum of the data (the model hrf+derivative does), and one can correct for this. An important aspect to consider beside the amplitude and time shift is the smoothness of images, as this impact statistics as well. Boosting relies on the ratio of parameters leading to rougher images. The recommendation is thus to smooth a little your data before running the GLM (e.g. 2mm) to increase SNR and then to re-estimate the response magnitude and smooth again (e.g. 6mm). At the group level, the smoothness achieved should match what you would have wanted if no re-estimation was performed.

 _spmup_hrf_boost:_

 This function re-estimate the response amplitude given the computed time shift, the time delay image is also returned. WARNING: the re-estimation is actually computed also if the 2nd derivative is present, and this is not recommended, unless a single condition is estimated, which means that by default the GLM should only include the hrf and 1st derivatives.

 _spmup_smooth_boostedfiles:_

 This function simply applies spm_smooth but re-masking before writing down the data. This is necessary as spmup_hrf_boost i) make images rougher and ii) introduces edges effect. The mask is simply taken from the computed GLM model.

## plots

_spmup_plotmotion:_

Plots separately translation, rotation and total displacement.

_spmup_plot_tsdiff:_

Plots the smalest and largest difference between time series in the time or frequency domain.

_MakeContours:_

from the current imaged data in spm_graphics, it makes an image of the contours of each blobs - this is useful to then show the raw statistical map and the signiticant areas highlighted by the contours.

_gp_parametric_plot:_

from a group level analysis on a parametric regressor, the average reponse (mean beta/con) is plotted along with the average model, if the model is of the same size between subjects.  All is returned in the workspace (with boostrapped confidence intervals).

_gp_event_plot:_

from a group level analysis, the the average reponse (mean beta/con) is plotted along with the average model, ie the reconstructed average hemodynamic reponse (also accounting for derivatives). All is returned in the workspace (with boostrapped confidence intervals).

## utilities

_spmup_auto_reorient:_

set the [0,0,0] coordinate automatically from one image to all images in a cell array - i.e. input the T1 and all fMRI data, the [0,0,0] coordinate is set for the T1 and this is applied to all fMRI data.

_spmup_auto_mask:_

Compute a mask from a time-series of memory mapped images, giving a similar (but more inclusive) mask than SPM.

_spmup_autocorrelation:_

This function computes (efficiently) the autocorrelation of all in mask  voxels of time series data.

_spmup_basics:_

Returns the most basic metrics: mean and std images from a time series.

_spmup_comp_dist2surf:_

compute the average distance from [0 0 0] to the surface of the brain (assuming at least segmented c1 c2 present, or better a surface gifti file).

_spmup_despike:_

'Despike' fMRI time series, by default using a median filter with a derived from the autocorrelation. Alternatively, user can input any of the matlab smooth function arguments.

_spmup_despike_reviewlog:_

QA for spmup_despike - indicating which voxels and how much despiking has been done.

_spmup_ernst:_

This returns the theoretical optimal flip angle for a given TR, i.e. cos(flip_angle) = exp(-TR/T1 time).

_spm_psc:_

Computes the percentage signal change scaled for the design matrix.

_spmup_resize:_

Routine to resize an image to match the desired bounding box and voxel size.

_spmup_roiextraction:_

Routine to extract the 1st eigen value for an arbitrary number of ROI and from an arbitrary number of images (like spm_voi).

_spmup_skip:_

Reads a 4D fmri data set and return a 4D set minus N initial volumes.
