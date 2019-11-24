function jobs = spmup_second_level_qa(varargin)

% routine calling spmup_find_dropout and find_group_outliers

spmup_find_dropout(varargin)

[bad_subject] = find_group_outliers(SPM_file)

