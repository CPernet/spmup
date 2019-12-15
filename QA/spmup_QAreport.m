function out = spmup_QAreport(varargin)

% meta function taking some of the outputs from subject level QA metrics
% to make a group level report and potentially flag outlying subjects
%
% FORMAT spmup_QAreport('metric',cellarray,'figure',fig_value)
%
% INPUTS - 'metric' can be 'anatQA', 'tSNR', 'volumecorr', 'spatialcorr' or 'displacement'
%        - cellarray is a cell array of the corresponding metrics, i.e.
%        a cell array of the anatQA or tSNR structures, a cell array of
%        volume correlation time courses, a cell array of slice_outliers, or
%        acell array of motion matrices ... note that for such arrays,
%        rows at taken to be subjects and columns to be sessions and/or runs
%        - 'figure' set to on (default) or save 
%
% OUTPUTS - a figure report 
%         - a structure with the same 'metric' input fields and corresponding
%         group level statistics
%
% Cyril Pernet 
% --------------------------------------------------------------------------
% Copyright (C) spmup team 2019

%% check inputs

out          = [];
fig          = 'on';
tSNR         = [];
anatQA       = [];
volumecorr   = [];
spatialcorr  = [];
displacement = [];

for in = 1:nargin
    if strcmpi(varargin{in},'anatQA')
        anatQA = varargin{in+1};
    elseif strcmpi(varargin{in},'volumecorr')
        volumecorr = varargin{in+1};
    elseif strcmpi(varargin{in},'spatialcorr')
        spatialcorr = varargin{in+1};
    elseif strcmpi(varargin{in},'displacement')
        displacement = varargin{in+1};
    elseif strcmpi(varargin{in},'tSNR')
        tSNR = varargin{in+1};
    elseif strcmpi(varargin{in},'figure')
        fig = varargin{in+1};
    else
        if ischar(varargin{in})
            warning('unknown argument %s, key ignored',varargin{in});
        end
    end
end

%% AnatQA

if ~isempty(anatQA)
    [n_subjects,p_repetitions] = size(anatQA);
    if n_subjects == 1
        anatQA = anatQA';
        [n_subjects,p_repetitions] = size(anatQA);
    end
    
    % make a matrix, get metrics and outliers, plot
    Fields = {'SNR','CNR','FBER','EFC','Asymmetry'};
    M = NaN(n_subjects,5,p_repetitions);
    for s = 1:n_subjects
        M(s,1,:) = anatQA{s}.SNR;
        M(s,2,:) = anatQA{s}.CNR;
        M(s,3,:) = anatQA{s}.FBER;
        M(s,4,:) = anatQA{s}.EFC;
        M(s,5,:) = anatQA{s}.asym;
    end
    
    fields_to_keep = sum(isnan(M))~=n_subjects;
    [Med,CI]       = HD(M(:,fields_to_keep)); 
    outliers       = spmup_comp_robust_outliers(M(:,fields_to_keep));
    
    if fields_to_keep(1) == 1
        out.anatQA.SNR.median   = Med(1);
        out.anatQA.SNR.CI       = CI(:,1)';
        out.anatQA.SNR.outliers = outliers(:,1)';
    end
    
    if fields_to_keep(2) == 1
        out.anatQA.CNR.median   = Med(2);
        out.anatQA.CNR.CI       = CI(:,2)';
        out.anatQA.CNR.outliers = outliers(:,2)';
    end
    
    if fields_to_keep(3) == 1
        out.anatQA.FBER.median   = Med(3);
        out.anatQA.FBER.CI       = CI(:,3)';
        out.anatQA.FBER.outliers = outliers(:,3)';
    end
    
    if fields_to_keep(4) == 1
        out.anatQA.EFC.median   = Med(4);
        out.anatQA.EFC.CI       = CI(:,4)';
        out.anatQA.EFC.outliers = outliers(:,4)';
    end
    
    if fields_to_keep(5) == 1
        out.anatQA.asym.median   = Med(5);
        out.anatQA.asym.CI       = CI(:,5)';
        out.anatQA.asym.outliers = outliers(:,5)';
    end

    make_figure('Summary Anat stats QA',fig,M,Fields,outliers)
end

%% tSNR

if ~isempty(tSNR)

end

%% volumecorr

if ~isempty(volumecorr)
    [n_subjects,p_repetitions] = size(volumecorr);
    if n_subjects == 1
        volumecorr = volumecorr';
        [n_subjects,p_repetitions] = size(volumecorr);
    end
    
   
    M = NaN(n_subjects,p_repetitions); % get median correlation
    D = NaN(n_subjects,p_repetitions); % get dispersion MAD   
    for p=1:p_repetitions
        for s = 1:n_subjects
            M(s,1:n(s,p)) = volumecorr{s,p_repetitions}';
        end
        
        [Med,CI]       = HD(M);
        outliers       = spmup_comp_robust_outliers(M);
        
         make_figure('Summary Anat stats QA',fig,M,Fields,outliers)
    end
    

end

%%
if ~isempty(spatialcorr)
    
end

%%
if ~isempty(displacement)
    
end

end

%% subfunction for median

function [HDQ,CIQ] = HD(X)

% Compute the Harrell-Davis estimate of the qth decile
% The vector x contains the data, and the desired decile is q
%
% FRANK E. HARRELL and C. E. DAVIS (1982).
% A new distribution-free quantile estimator
% Biometrika 69 (3): 635-640. doi: 10.1093/biomet/69.3.635

[p,N]  = size(X); % number of estimates to compute
q      = .5; % median
nboot  = 1000; % use for 95% CI
alphav = 5/100; % for percentile bootstrap when n<11
rng shuffle

%% compute
table = randi(p,p,nboot);
for i=1:N
    nanindex = isnan(X(:,i));
    x        = X(~nanindex,i);
    HDQ(i)   = get_HD(x,q);
    
    % The constant c was determined so that the
    % probability coverage of the confidence interval is
    % approximately 95% when sampling from normal and
    % non-normal distributions
    n=length(x);
    
    if n<=10 % hd estimate of the deciles cannot be computed, using percentile boostrap
        for kk=1:nboot
            if sum(nanindex) ~=0
                values = find(nanindex);
                tmp = table(:,kk);
                for l=1:length(values)
                    tmp(tmp==values(l))=[];
                end
                D = X(tmp,i);
            else
                D = X(table(:,kk),i) ; % applies the sample resampling for each column
            end
            Mb(kk) = HD(D,.5);
        end
        
        Mb = sort(Mb);
        Mb(isnan(Mb)) = [];
        low = round((alphav*length(Mb))/2);
        high = length(Mb) - low;
        CIQ(1,i) = Mb(low+1);
        CIQ(2,i) = Mb(high);
        
    else
        if n <=21 &&  q<=.2 || n <=21 && q>=.8
            c = -6.23./n+5.01;
        elseif n<=40 && q<=.1 || n<=40 && q>=.9
            c = 36.2./n+1.31;
        else
            c = 1.96 + .5064.* (n.^-.25);
        end
        
        % do bootstrap
        for kk=1:nboot
            boot(kk)=get_HD(randsample(x,n,true),q);
        end
        bse = std(boot,0); % normalize by (n-1)
        CIQ(1,i) = HDQ(i)-c.*bse;
        CIQ(2,i) = HDQ(i)+c.*bse;
    end
    clear x n c boot bse
end

end

function HD = get_HD(x,q)
% that's the Harrell-Davis estimator
n   = length(x);
m1  = (n+1).*q;
m2  = (n+1).*(1-q);
vec = 1:n;
w   = betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2);
y   = sort(x);
HD  = sum(w(:).*y(:));
clear vec w y
end

function make_figure(fig_name,fig,M,Fields,outliers)

figure('Name',fig_name)
if strcmpi(fig,'on')
    set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
else
    set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1],'visible','off')
end

x = (rand(size(M,1),1).*0.3)+0.85;
for m=1:5
    subplot(1,5,m); boxplot(M(:,m),'boxstyle','outline', 'whisker',1.5,'widths',0.5);
    a = findall(gca,'type','line');
    for l=1:size(a,1)
        if strcmp(a(l).LineStyle,'-') || strcmp(a(l).LineStyle,'--')
            set(a(l), 'linewidth',2); set(a(l), 'Color', [0 0 0]);
        end
    end
    set(a(2), 'Color', [1 0 0]);
    title(Fields{m}); box on; grid on;
    h = findobj(gca,'Tag','Box');
    for j=length(h):-1:1
        patch(get(h(j),'XData'),get(h(j),'YData'),[0 0 0],'FaceAlpha',.4); % fill with transparent color
    end
    hold on; plot(x,M(:,m), 'ko','linewidth', 1);
    if m <=size(outliers,2)
        tmp = M(:,m).*outliers(:,m); tmp(tmp==0)=NaN;
        plot(x,tmp,'ro','linewidth', 2); clear tmp
    end
end

if strcmpi(fig,'save')
    print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(pwd,'spmup_QAreport.ps'));
    close(gcf)
end
end

