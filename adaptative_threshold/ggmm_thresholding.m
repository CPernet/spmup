function [threshold,bic] = ggmm_thresholding(stat_filename, mask_filename)

% read the data (spmT and mask), call the EM code and return the T value
% for the cluster forming threshold
%
% FORMAT: [threshold,em,bic] = ggmm_thresholding(stat_filename, mask_filename)
%
% INPUT: stat_filename: full name (ie with path) of the spmT image
%        mask_filename: full name (ie with path) of the mask image
% 
% OUTPUT: threshold: NaN if Gaussian only was found best model
%                    one value X if positive Gamma was found best model,
%                    two values X1 X2 if negative and positive Gamma was found best model,
%         bic: Bayesian information criteria of each model
%
% Coded by Chris Gorgolewski 29 June 2012
% Edited (help, outputs, etc) Cyril Pernet 4 July 2012
% ------------------------------------------
% This file is part of the adaptative_threshold SPM plug-in
% Copyright (C) Cyril Pernet and Chris Gorgolewski 2012


    mask_data = spm_read_vols(spm_vol(mask_filename));
    stat_data = spm_read_vols(spm_vol(stat_filename));
    stat_data = stat_data(mask_data > 0);
    no_signal_components = {Gaussian(0, 10)};
    noise_and_activation_components = {Gaussian(0, 10), Gamma(4, 5, 0)};
    noise_activation_and_deactivation_components = {Gaussian(0, 10), Gamma(4, 5, 0), NegativeGamma(4,5, 0)};

    models = {no_signal_components, noise_and_activation_components, noise_activation_and_deactivation_components};
    best_bic = inf;

    figure('Name','Gamma-Gaussian mixture model of T voxel values','Visible','On');
    set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized', 'outerposition',[0 0 1 1]); 
    for i=1:length(models)
        em = ExpectationMaximization(models{i});
        em.fit(stat_data);
        subplot(1,3,i);
        em.plot(stat_data);
        drawnow
        bic(i) = em.BIC(stat_data);
        if bic(i) < best_bic
            best_bic = bic(i);
            best_model = em;
        end
    end

    if size(best_model.components,2) == 1
        fprintf('Best model is Gaussian only, ie no signal was found!\n')
        threshold = NaN; 
    elseif size(best_model.components,2) == 2
        fprintf('Best model is Gaussian + Positive Gamma\n')
        pp = best_model.posteriors(stat_data);
        active_map = (pp(:, 2) > pp(:, 1)) & stat_data > 0.01;
        threshold = min(stat_data(active_map));
    elseif size(best_model.components,2) == 3
        fprintf('Best model is Negative Gamma + Gaussian + Positive Gamma\n')
        pp = best_model.posteriors(stat_data);
        active_map = (pp(:, 2) > pp(:, 1)) & stat_data > 0.01;
        deactive_map = (pp(:, 3) > pp(:, 1)) & stat_data < 0.01;
        threshold = [max(stat_data(deactive_map)) min(stat_data(active_map))];
    end

end