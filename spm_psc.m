function [SF]=spm_psc(SPMs,scaling_factor)

% either return the maximum heights of regressors among SPMs or
% computes the percentage signal change using the information from the
% SPM.mat - if the scaling factor is not specifed, figures will pop out for
% each subject and the SF will have to be specified.
% 
% FORMAT: [heights] = spm_psc 
%         [heights] = spm_psc(SPMs)
%         spm_psc(SPM,scaling_factor)
% 
% INPUT: SPMs is a cell array or structure of SPM.mat files
%        scaling factor is the scaling to apply to compute PSC
%
% OUTPUT: [heights] returns the maximum values for each condition in each
%         subject - this max can be used as scaling factor
%         none = write down images of PSC
%
% Cyril Pernet
% -------------

%% check varargin
SF = [];
if nargin == 0 % pick up as many SPM.mat files as needed
    [SPMs,sts] = spm_select(Inf,'mat','Select SPM.mat files');
    if sts == 0
        return
    end
    
elseif nargin == 2
    SF = varargin{2};
end


% --------------------------
% loop per subject

for subject = 1:size(SPMs,1)
    fprintf('processing subject %g \n',subject)
    try
        [path,file,ext]= fileparts(SPMs{subject});
    catch SPMs_file_issue
        [path,file,ext]= fileparts(SPMs(subject,:));
    end
    cd(path); load(file)
    
    %% 1 get the max of ideal trial in the super sampled design matrix
    xBF.order = SPM.xBF.order;
    if isempty(SF)
        xBF.dt = SPM.xBF.dt;
        xBF.name = SPM.xBF.name;
        xBF.length = SPM.xBF.length;
        xBF = spm_get_bf(xBF);
        if strcmp(SPM.xBF.name,'hrf (with time derivative)') 
            event = xBF.bf*[1 1]';
        elseif strcmp(SPM.xBF.name,'hrf (with time and dispersion derivatives)')
            event = xBF.bf*[1 1 1]';
        else
            event = xBF.bf;
        end
        SF = max(event);
    end
    
    %% 2 if derivative(s) present - compute the boosted hrf
    if xBF.order>1 && ~isdir('hrf_boost')
        spmup_hrf_boost([SPM.swd filesep 'SPM.mat'])
    end
    
    fprintf('computing and writing PSC images using a value of %g as scaling factor \n',SF)
    
    %% 3 compute PSC using the hrf or boosted hrf
    mkdir PSC
    s = size(SPM.Sess,2);
    for session=1:s
        U = spm_get_ons(SPM,s);
        [X,Xn,Fc] = spm_Volterra(U,xBF.bf,SPM.xBF.Volterra); % re-create the super-sampled design
        nb_conditions = size(SPM.Sess(s).U,2);
        for c=1:nb_conditions
            columns = SPM.Sess.Fc(c).i; % which columns for this condition
            hrf_param = [1:3:length(columns)]; % in case we have parametric regressors
            for h=1:length(hrf_param)
                regressors = columns(hrf_param(h) + (0:(xBF.order-1)));
                combined_regressor = X(:,regressors)*ones(xBF.order,1); % using all functions this is the hrf model

                for l=1:xBF.order % load beta files
                    if regressors(l)<10
                        name = ['beta_000' num2str(regressors(l)) '.img,1'];
                    else
                        name = ['beta_00' num2str(regressors(l)) '.img,1'];
                    end
                    beta{l} = spm_read_vols(spm_vol([pwd filesep name]));
                end
                % now combine beta values - scaled by the new regressor
%                 H  = sqrt(((beta{1}.*beta{1}).*sum(X(:,regressors(1)).^2))+...
%                     ((beta{2}.*beta{2}).*sum(X(:,regressors(2)).^2))+...
%                     ((beta{3}.*beta{3}).*sum(X(:,regressors(3)).^2)));
                H = 0;
                for l=1:xBF.order
                    H = H + ((beta{l}.*beta{1}).*sum(X(:,regressors(l)).^2));
                end
                H = sqrt(H);
                % keep the sign of beta hrf
                beta_sign = beta{1} ./ abs(beta{1});
                
                % finally get the PSC
                last_column = size(SPM.xX.X,2);
                if last_column<10
                    name = ['beta_000' num2str(last_column) '.img,1'];
                else
                    name = ['beta_00' num2str(last_column) '.img,1'];
                end
                V = spm_vol([pwd filesep name]);
                constant = spm_read_vols(V);
                PSC = beta_sign .* (H.*SF./constant.*100);
                % write as a matrix and as an image using column number
                if regressors(1)< 10
                    V.fname = [pwd filesep 'PSC' filesep 'PSC_000' num2str(regressors(1)) '.img'];
                else
                    V.fname = [pwd filesep 'PSC' filesep 'PSC_00'  num2str(regressors(1)) '.img'];
                end
                V.descrip = 'Percentage signal change image';
                V.private.descrip = 'Percentage signal change image';
                spm_write_vol(V,PSC);
            end
        end
    end
end

