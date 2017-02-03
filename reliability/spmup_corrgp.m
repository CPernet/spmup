function [W,B,Stats] = spmup_corrgp(series1,series2,masks,figout,threshold)

% Routine that allows to test for image reliability based on correlations. 
% The function calls iteratively spmup_corr to compute the Pearson and
% concordance correlations between images from series1 and series2
% All pairs are computed to produce within vs. between pairs estimates
% For instance for series1 = [A1 A2 A3] and series2 = [B1 B2 B3],
% Pearson and concordance correlations are computed for [A1B1] [A2B2] and
% [A3B3] giving the within pairs correlations, and mean([A1B2 A1B3]),
% mean([A2B1 A2B3]) and mean([A3B1] [A3B2]) giving the between pairs
% correlations. A bootstrap is then computed to test if the diff = 0.
%
% FORMAT = [W,B,Stats] = spmup_corrgp(series1,series2,masks,figout,threshold)
%
% INPUT series1 series of image filename (see spm_select)
%       series2 series of image filename (see spm_select)
%       masks series of image filename (see spm_select)
%       figout 1/0 (default) to get all correlation figures out
%       threshold (optional) if masks are not binary, threshold to apply
%
% output W are the within pairs correlations (Pearson and Concordance)
%        B are the within pairs correlations (Pearson and Concordance)
%        Stats is structure
%
% Cyril Pernet
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2017

if nargin == 3
    figout = 0;
    threshold = 0;
end

if size(series1,1) ~= size(series2,1)
    error('series 1 and 2 are of different size')
end

if size(masks,1) == 1
    disp('only one mask used for all images')
    masks = repmat(mask,size(series1,1),1);
end

%% do the within loop
N = size(series1,1); W = nan(N,2); 
for n=1:N
    fprintf('Computing within pairs correlations: pair %g/%g \n',n,N)
    [W(n,1),W(n,2)]=spmup_corr(series1(n,:),series2(n,:),masks(n,:),figout,threshold);
end
    
%% do the between loop
combinations = nchoosek([1:size(series1,1)],2);
MP = NaN(N,N); MC = NaN(N,N); 
for n=1:length(combinations)
    fprintf('Computing between pairs correlations: pair %g/%g \n',n,length(combinations))    
    [MP(combinations(n,1),combinations(n,2)),MC(combinations(n,1),combinations(n,2))] = spmup_corr(...
        series1(combinations(n,1),:),series2(combinations(n,2),:),masks(n,:),figout,threshold);
    MP(combinations(n,2),combinations(n,1)) = MP(combinations(n,1),combinations(n,2));
    MC(combinations(n,2),combinations(n,1)) = MC(combinations(n,1),combinations(n,2));
end
B = [nanmean(MP,2) nanmean(MC,2)];


%% do the stats
if nargout == 3
   Stats.means_within    = mean(W);
   Stats.means_between   = mean(B);
   Stats.medians_within  = median(W);
   Stats.medians_between = median(B);

   % percentile bootstrap on the difference
   nboot = 600; alphav = 5/100;
   low = round((alphav*nboot)/2); high = nboot - low;
   WP = W(:,1); WP = WP(randi(size(W,1),size(W,1),nboot));
   BP = B(:,1); BP = BP(randi(size(B,1),size(B,1),nboot));
   WC = W(:,2); WC = WC(randi(size(W,2),size(W,2),nboot));
   BC = B(:,2); BC = BC(randi(size(B,2),size(B,2),nboot));
   
   % mean 
   MWP = sort(mean(WP,1)); MWC = sort(mean(WC,1)); 
   Stats.CI.means_within(1,:) = [MWP(low)  MWC(low)];
   Stats.CI.means_within(2,:) = [MWP(high) MWC(high)];
   clear MWP MWC
   
   MBP = sort(mean(BP,1)); MBC = sort(mean(BC,1)); 
   Stats.CI.means_between(1,:) = [MBP(low)  MBP(low)];
   Stats.CI.means_between(2,:) = [MBP(high) MBP(high)];
   clear MBP MBC
   
   DP = sort(mean(WP,1)-mean(BP,1)); 
   DC = sort(mean(WC,1)-mean(BC,1)); 
   Stats.means_difference = mean(W)-mean(B);
   Stats.CI.means_difference(1,:) = [DP(low)  DC(low)];
   Stats.CI.means_difference(2,:) = [DP(high) DC(high)];
   pb = sum(DP > 0) / nboot; 
   Stats.pval.means_difference = 2*min(pb,1-pb);
   pb = sum(DC > 0) / nboot; 
   Stats.pval.means_difference(2) = 2*min(pb,1-pb);
   clear DP DC
  
   % median 
   MWP = sort(median(WP,1)); MWC = sort(median(WC,1)); 
   Stats.CI.medians_within(1,:) = [MWP(low)  MWC(low)];
   Stats.CI.medians_within(2,:) = [MWP(high) MWC(high)];
   clear MWP MWC
   
   MBP = sort(median(BP,1)); MBC = sort(median(BC,1)); 
   Stats.CI.medians_between(1,:) = [MBP(low)  MBP(low)];
   Stats.CI.medians_between(2,:) = [MBP(high) MBP(high)];
   clear MBP MBC
   
   DP = sort(median(WP,1)-median(BP,1)); 
   DC = sort(median(WC,1)-median(BC,1)); 
   Stats.medians_difference = median(W)-median(B);
   Stats.CI.medians_difference(1,:) = [DP(low)  DC(low)];
   Stats.CI.medians_difference(2,:) = [DP(high) DC(high)];
   pb = sum(DP > 0) / nboot; 
   Stats.pval.medians_difference = 2*min(pb,1-pb);
   pb = sum(DC > 0) / nboot; 
   Stats.pval.medians_difference(2) = 2*min(pb,1-pb);
   clear DP DC   
   
end





