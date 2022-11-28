function status = spmup_check_matlab_version()
    status = 0;
    % matlab 2019a
    % see https://en.wikipedia.org/wiki/MATLAB#Release_history
    matlabVersion = '9.6.0'; 
    if verLessThan('matlab', matlabVersion)
        warning('Some features are not implemented for MATLAB version < R2020a.');
        status = 1;
    end
end
