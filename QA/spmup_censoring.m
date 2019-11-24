function design = spmup_censoring(varargin)

% routine that computes robust outliers for each column of the data in
% and return a matrix of censsoring regressors (0s and a 1 per column)
%
% FORMAT censoring_regressors = spmup_censoring(realignemt.txt,data)
%        censoring_regressors = spmup_censoring(realignemt.txt,data,'Voltera','on/off')
%
% INPUT  realignment.txt is a motion parameter file 
%        data is a n volumes * m matrix to censor column wise
%        Voltera indicate to expand the motion parameters
%
% OUTPUT censoring_regressors a matrix with ones in each column for
%        outliers found in coumns of data (also saved as design.txt)
% 
% Exemple of data matrix data = [outlying_voxels' r_course']
% Volterra expansion of motion parameters 
% Friston et al (1996) Magn Reson Med. 1996;35:346:355.

% Cyril Pernet
% --------------------------------
% Copyright (C) spmup 2017

%% inputs

Voltera = 'off';
if nargin >2
    for v=3:nargin
        if strcmpi(varargin{v},'Voltera')
            Voltera = varargin{v+1};
        end
    end
end
data = varargin{2};

% motion

[filepath,filename]=fileparts(varargin{1});
motion = load(varargin{1}); % load data
   
if strcmpi(Voltera,'on')
    D      = diff(motion,1,1); % 1st order derivative
    D      = [zeros(1,6); D];  % set first row to 0
    motion = [motion D motion.^2 D.^2];
end

%% get outliers

if ~isempty(data)
    k = sqrt(chi2inv(0.975,1));
    [n,p] = size(data);
    distance = NaN(n,p);
    for p=1:size(data,2)
        tmp = data(:,p);
        points = find(~isnan(tmp));
        tmp(isnan(tmp)) = [];
        
        % compute all distances
        n = length(tmp);
        for i=1:n
            j = points(i);
            indices = [1:n]; indices(i) = [];
            distance(j,p) = median(abs(tmp(i) - tmp(indices)));
        end
        
        % get the S estimator
        % consistency factor c = 1.1926;
        Sn = 1.1926*median(distance(points,p));
        
        % get the outliers in a normal distribution
        outliers(:,p) = (distance(:,p) ./ Sn) > k; % no scaling needed as S estimates already std(data)
        outliers(:,p) = outliers(:,p)+isnan(data(:,p));
    end
    
    % create the motion censoring regressors simply diagonalizing each outlier
    % columns and then removing empty columns
    censoring_regressors = [];
    for column = p:-1:1
        censoring_regressors = [censoring_regressors diag(outliers(:,column))];
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
    design = [motion censoring_regressors];
else
    design = motion;
end
save(fullfile(filepath,[filename(4:end) '_design.txt']),'design','-ascii')

% % quick check
% if sum(censoring_regressors(:))/size(censoring_regressors,1)*100 > 10
%     warndlg('Censoring was performed but more than 10% of data were marked which can be problematic fitting data');
% end

