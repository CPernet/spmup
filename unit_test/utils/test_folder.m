function value = test_folder()
 value = spm_file(fullfile(fileparts(mfilename('fullpath')), '..'), 'cpath');
end