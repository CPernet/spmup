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
% Ref.  Pernet CR (2014) Misconceptions in the use of the General Linear
% Model applied to functional MRI: a tutorial for junior neuro-imagers.
% Front. Neurosci. 8:1. doi: 10.3389/fnins.2014.00001
% <https://www.frontiersin.org/articles/10.3389/fnins.2014.00001/full>
%
% Cyril Pernet - University of Edinburgh
% Thanks to Robin Ince for improving the code
% -----------------------------------------
% Copyright (c) SPM Utility Plus toolbox

%% check varargin
if nargin == 0 % pick up as many SPM.mat files as needed
    [SPMs,sts] = spm_select(Inf,'mat','Select SPM.mat files');
    if sts == 0
        return
    end
    
    % get SF
    option = questdlg('scale for a typical trial or for this design?','option','trial','design','trial');
    SF     = spmup_psc(SPMs,option);
    if length(SF) == 1
        SF = repmat(SF,[1,size(SPMs,1)]);
    end
    clear option
    
else
    if nargout == 1
        option = varargin{1};
    else
        SF = varargin{1};
        if length(SF) == 1
            SF = repmat(SF,[1,size(SPMs,1)]);
        end
    end
end

%% loop per subject
for subject = 1:size(SPMs,1)
    
    % read the design matrix
    if iscell(SPMs)
        [path,file,ext]= fileparts(SPMs{subject});
    else
        [path,file,ext]= fileparts(SPMs(subject,:));
    end
    SPM        = load(fullfile(path,[file ext]));
    SPM        = SPM.SPM;
    xBF.dt     = SPM.xBF.dt;
    xBF.name   = SPM.xBF.name;
    xBF.length = SPM.xBF.length;
    xBF.order  = SPM.xBF.order;
    xBF        = spm_get_bf(xBF); % rebuild the hrf model used
    
    %% Get the Scaling Factor
    if exist('option','var')
        
        if strcmpi(option,'trial')
            if strcmp(SPM.xBF.name,'hrf (with time derivative)')
                event = xBF.bf*[1 1]';
            elseif strcmp(SPM.xBF.name,'hrf (with time and dispersion derivatives)')
                event = xBF.bf*[1 1 1]';
            else
                event = xBF.bf;
            end
            scaling_factor(subject) = max(event);
            
        elseif strcmpi(option,'design')
            s = size(SPM.Sess,2);
            for session= s:-1:1
                U = spm_get_ons(SPM,s);
                X = spm_Volterra(U,xBF.bf,SPM.xBF.Volterra); % re-create the super-sampled design
                if strcmpi(option,'design')
                    sf_sess(session) = max(X(:));
                end
            end
            scaling_factor(subject) = max(sf_sess);
            
        else
            error('unkown option for scaling factor')
        end
    end
    
    %% 2 compute PSC 
    
    if exist('SF','var')
        fprintf('computing and writing PSC images using a value of %g as scaling factor \n',SF(subject));
        
        %% 3 compute PSC using the hrf or boosted hrf
        mkdir(fullfile(path,'PSC'))
        s = size(SPM.Sess,2);
        hrf_indices = [];
        for session= 1:s
            U             = spm_get_ons(SPM,s);
            X             = spm_Volterra(U,xBF.bf,SPM.xBF.Volterra); % re-create the super-sampled design
            nb_conditions = size(SPM.Sess(s).U,2);
            for c=1:nb_conditions
                columns   = SPM.Sess(session).Fc(c).i; % which columns for this condition
                hrf_param = 1:3:length(columns); % in case we have parametric regressors
                for h=1:length(hrf_param)
                    hrf_indices = [hrf_indices SPM.Sess(session).col(columns(hrf_param(h)))]; % keep for contrasts
                    regressors  = columns(hrf_param(h) + (0:(xBF.order-1)));
                    for l=1:xBF.order % load beta files
                        name    = [SPM.Vbeta(regressors(l)).fname ',1'];
                        beta{l} = spm_read_vols(spm_vol([SPM.swd filesep name]));
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
                    name      = [SPM.Vbeta(max(SPM.Sess(s).col)+session).fname]; % constant
                    [~,~,ext] = fileparts(name); % nii or img
                    V         = spm_vol([SPM.swd filesep name]);
                    constant  = spm_read_vols(V);
                    PSC       = beta_sign .* (H.*SF(subject)./constant.*100);
                    % write as a matrix and as an image using column number
                    V.fname           = fullfile(SPM.swd,'PSC',[sprintf('PSC_%04d',SPM.Sess(session).col(columns(hrf_param(h)))) ext]);
                    V.descrip         = 'Percentage signal change image';
                    V.private.descrip = 'Percentage signal change image';
                    spm_write_vol(V,PSC);
                end
            end
        end
    end
end

% cleanup SF
if exist('scaling_factor','var')
    if length(unique(scaling_factor)) == 1
        SF = scaling_factor(1);
    else
        warning('different scaling factors obtained for each design')
        SF = scaling_factor;
    end
end

