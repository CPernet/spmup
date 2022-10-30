function spmup_setparallel(varargin)

% simple routine to set the up the parallel processing
%
% FORMAT spmup_setparallel(N)
% set N if you want to use a specific number of processors 
% useful on servers with dynamic attribution to limit oneself
%
% Cyril Pernet & Remi Gau
% --------------------------
%  Copyright (C) SPMUP Team 

if nargin == 0
    N = [];
else
    N = varargin{1};
end

if evalin( 'base', 'exist(''options'',''var'') == 1' )
    options = evalin('base','options');
    N       = options.Ncores;
end

addons = ver;
if any(any(contains(struct2cell(addons), 'Parallel Computing Toolbox')))
    p = gcp('nocreate');
    if isempty(p) % i.e. the parallel toolbox is not already on
        if isempty(N)
            
            % check how many cores to use
            % ---------------------------
            N = getenv('NUMBER_OF_PROCESSORS'); % logical number of cores (i.e. count hyperthreading)
            if isempty(N)
                N = feature('numcores');          % physical number of cores
            end
            
            if ischar(N)
                N = str2double(N);
            end
            
            % check and set the local profile NumWorkers
            % ----------------------------------
            c            = parcluster;
            c.NumWorkers = N;
            saveProfile(c);
            
            % go
            % --
            try
                parpool(N-1);
            catch flasestart %#ok<NASGU> 
                delete(c.Jobs)
                parcluster('local')
                parpool(N-1);
            end
        else
            try
                parpool(N);
            catch
                disp(['Parallel computing could not be set up.', ...
                      ' Nothing to worry about (except slower computation in some cases)']);
            end
        end
    end
else
    disp(['Parallel computing could not be set up.', ...
          ' Nothing to worry about (except slower computation in some cases)']);
end
