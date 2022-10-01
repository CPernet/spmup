function spmup_setparallel

% simple routine to set the up the parallel processing
%
% Cyril Pernet & Remi Gau
% --------------------------
%  Copyright (C) SPMUP Team 


N      = 10;
addons = ver;
if any(strcmpi('Parallel Computing Toolbox',arrayfun(@(x) x.Name, addons, "UniformOutput",false)))
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
            % saveProfile(c);
            
            % go
            % --
            parpool(N-1);
        else
            parpool(N);
        end
    end
else
    disp('Parallel toolbox not found - nothing to worry about (except slower computation in some cases)');
end
