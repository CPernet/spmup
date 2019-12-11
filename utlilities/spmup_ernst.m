function flip_angle = spmup_ernst(TR,ToneTime)

% this returns the theoretical optimal flip angle for a given TR:
% cos(flip_angle) = exp(-TR/T1)
%
% FORMAT flip_angle = spmup_ernst(TR,ToneTime)
%
% INPUT TR of the pulse sequence
%       ToneTime (optional) is the T1 recovery time, default is 1400 optimized for gray matter

% Cyril Pernet August 2017
% --------------------------------------------------------------------------
% Copyright (c) SPM Utility Plus toolbox

if nargin == 1
  ToneTime = 1400;
end

flip_angle = acosd(exp(-TR/ToneTime));




