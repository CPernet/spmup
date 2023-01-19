function [IS_GITHUB, pth] = is_github_ci(verbose)
  %

  % (C) Copyright 2021 Remi Gau

  if nargin < 1 || isempty(verbose)

    verbose = false;

  end

  IS_GITHUB = false;

  GITHUB_WORKSPACE = getenv('HOME');

  if strcmp(GITHUB_WORKSPACE, '/home/runner')

    if verbose
      fprintf(1, '\n WE ARE RUNNING IN GITHUB CI\n');
    end

    IS_GITHUB = true;
    pth = GITHUB_WORKSPACE;

  end

end
