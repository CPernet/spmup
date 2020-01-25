function SF = spm_psc(SPMs,varargin)

% either return the maximum heights of regressors among SPMs or
% computes the percentage signal change using the information from the
% SPM.mat - if the scaling factor is not specifed, figures will pop out for
% each subject and the SF will have to be specified.
%
% FORMAT: spm_psc
%         Scaling_Factor = spm_psc(SPMs,option)
%         spm_psc(SPM,scaling_factor)
%
% INPUT: if no inputs, do it all and call GUI
%        - SPMs is a cell array or structure of SPM.mat files
%        - option is either 'trial' to get the scaling factor for a typical
%        hrf (given your model parameters) or 'design' to get a scaling
%        for the design at hand (ie taking the max)
%        - scaling factor a single value or a vector of scaling value to
%          apply to compute PSC
%
% OUTPUT: Scaling Factor (SF) returns the value(s) to use for PSC
%         if no output specify, writes down images of PSC
%
% Reference: Pernet (2014) Misconceptions in the use of the GLM ...
% <https://www.frontiersin.org/articles/10.3389/fnins.2014.00001/full>
%
% Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

%% check varargin
if nargin == 0 % pick up as many SPM.mat files as needed
    [SPMs,sts] = spm_select(Inf,'mat','Select SPM.mat files');
    if sts == 0
        return
    end
    option = questdlg('scale for a typical trial or for this design?','option','trial','design','trial');
    SF = spm_psc(SPMs,option);
    if length(unique(SF)) > 1
        warning('different scaling factors obtained for each design')
    end
    
elseif nargin == 2
    if nargout == 1
        option = varargin{2};
    else
        SF = varargin{2};
    end
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
    SPM = load(fullfile(path,file));
    SPM = SPM.SPM;
    
    %% 1 get the max of ideal trial in the super sampled design matrix
    
    xBF.order = SPM.xBF.order;
    if strcmpi(option,'trial')
        xBF.dt     = SPM.xBF.dt;
        xBF.name   = SPM.xBF.name;
        xBF.length = SPM.xBF.length;
        xBF        = spm_get_bf(xBF); % rebuild the hrf model used
        if strcmp(SPM.xBF.name,'hrf (with time derivative)')
            event = xBF.bf*[1 1]';
        elseif strcmp(SPM.xBF.name,'hrf (with time and dispersion derivatives)')
            event = xBF.bf*[1 1 1]';
        else
            event = xBF.bf;
        end
        SF = max(event);
    elseif strcmpi(option,'design')
        s = size(SPM.Sess,2);
        for session= s:-1:1
            U = spm_get_ons(SPM,s);
            X = spm_Volterra(U,xBF.bf,SPM.xBF.Volterra); % re-create the super-sampled design
            if strcmpi(option,'design')
                sf_sess(session) = max(X(:));
            end
        end
        SF(subject) = max(sf_sess); % update
    else
        error('unkown option for scaling factor')
    end
    
    %% 2 if derivative(s) present - compute the boosted hrf
    if xBF.order>1 && ~isfolder('hrf_boost')
        spmup_hrf_boost([SPM.swd filesep 'SPM.mat'])
    end
    fprintf('computing and writing PSC images using a value of %g as scaling factor \n',SF);
    
    %% 3 compute PSC using the hrf or boosted hrf
    mkdir(fullfile(path,'PSC'))
    s = size(SPM.Sess,2);
    for session= 1:s
        U = spm_get_ons(SPM,s);
        X = spm_Volterra(U,xBF.bf,SPM.xBF.Volterra); % re-create the super-sampled design
        nb_conditions = size(SPM.Sess(s).U,2);
        for c=1:nb_conditions
            columns = SPM.Sess(session).Fc(c).i; % which columns for this condition
            hrf_param = [1:3:length(columns)]; % in case we have parametric regressors
            for h=1:length(hrf_param)
                regressors = columns(hrf_param(h) + (0:(xBF.order-1)));
                combined_regressor = X(:,regressors)*ones(xBF.order,1); % using all functions this is the hrf model                for l=1:3 % load beta files
                for l=1:xBF.order % load beta files
                    name = [SPM.Vbeta(regressors(l)).fname ',1'];
                    beta{l} = spm_read_vols(spm_vol([pwd filesep name]));
                end
                
                % now combine beta values - scaled by the new regressor
                H = 0;
                for l=1:xBF.order
                    H = H + ((beta{l}.*beta{l}).*sum(X(:,regressors(l)).^2));
                end
                H = sqrt(H);
                % keep the sign of beta hrf
                beta_sign = beta{1} ./ abs(beta{1});
                
                % finally get the PSC
                last_column = size(SPM.xX.X,2)-s+session; % 1 cst per session
                name = [SPM.Vbeta(last_column).fname ',1'];
                [~,~,ext] = fileparts(SPM.Vbeta(last_column).fname); % nii or img
                V        = spm_vol([pwd filesep name]);
                constant = spm_read_vols(V);
                PSC      = beta_sign .* (H.*SF./constant.*100);
                % write as a matrix and as an image using column number
                V.fname = fullfile(pwd,'PSC',[sprintf('PSC_%04d',regressors(1)) ext]);
                V.descrip         = 'Percentage signal change image';
                V.private.descrip = 'Percentage signal change image';
                spm_write_vol(V,PSC);
            end
        end
    end
end

%% combine PSC maps using contrast weights if any


if size(Scaling_Factor,2) == 1
    fprintf('PSC done - Scaling Factor %g: \n',Scaling_Factor)
else
    fprintf('PSC done - Scaling Factors: \n')
    for s=1:size(Scaling_Factor,2)
        fprintf('subject %g = %g \n',s,Scaling_Factor(s));
    end
end

