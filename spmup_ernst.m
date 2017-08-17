function flip_angle = spmup_ernst(TR)

% assuming T1 ~ 1400ms for gray matter, this returns the optimal flip
% angle for a given TR: cos(flip_angle) = exp(-TR/T1)
%
% Cyril Pernet August 2017
% --------------------------------------------------------------------------
% Copyright (c) SPM Utility Plus toolbox

flip_angle = acosd(exp(-TR/1400));




