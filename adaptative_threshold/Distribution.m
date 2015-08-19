% update the object with the type of distribution to fit
% ------------------------------------------
% This file is part of the adaptative_threshold SPM plug-in
% Copyright (C) Cyril Pernet and Chris Gorgolewski 2012

classdef Distribution < handle
    properties(Constant,Abstract)
        free_parameters
    end
    methods(Abstract)
        fit_weighted(obj, data, weights);
        pdf(obj, data);
    end 
end

