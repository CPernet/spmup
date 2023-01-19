function [design,censoring_regressors] = spmup_censoring(data)

% routine that computes robust outliers for each column of the data in
% and return a matrix of censsoring regressors (0s and a 1 per column)
%
% FORMAT [design,censoring_regressors] = spmup_censoring(data)
%        [design,censoring_regressors] = spmup_censoring(data)
%
% INPUT  data is a n volumes * m matrix to censor column wise
%
% OUTPUT design is the augmented design matrix of regressors
%        censoring_regressors matrix with ones in each column for
%        outliers found in coumns of data 
%        if no output specified saves it as design.txt
% 
% Exemple of data matrix data = [outlying_voxels' r_course']
% Friston et al (1996) Magn Reson Med. 1996;35:346:355.
%
% Cyril Pernet
% --------------------------
%  Copyright (C) SPMUP Team 


%% get outliers
design = [];
censoring_regressors = [];

if ~isempty(data)
    outliers = spmup_comp_robust_outliers(data,'Carling');

    % create the motion censoring regressors simply diagonalizing each outlier
    % columns and then removing empty columns
    censoring_regressors = [];
    for column = size(outliers,2):-1:1
        censoring_regressors = [censoring_regressors diag(outliers(:,column))]; %#ok<*AGROW>
    end
    censoring_regressors(:,sum(censoring_regressors,1)==0) = []; % remove empty columns
    % check a volume is not indexed twice
    for r=1:size(censoring_regressors,1)
        rep = find(censoring_regressors(r,:));
        if length(rep)>1
            censoring_regressors(r,rep(2:end)) = 0;
        end
    end
    censoring_regressors(:,sum(censoring_regressors,1)==0) = [];
    design = [data censoring_regressors];
end

if nargout == 0 && ~isempty(design)
    save(fullfile(filepath,[filename(4:end) '_design.txt']),'design','-ascii')
end

% % quick check
% if sum(censoring_regressors(:))/size(censoring_regressors,1)*100 > 10
%     warndlg('Censoring was performed but more than 10% of data were marked which can be problematic fitting data');
% end

