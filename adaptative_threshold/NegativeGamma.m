% create the negative Gamma model
% ------------------------------------------
% This file is part of the adaptative_threshold SPM plug-in
% Copyright (C) Cyril Pernet and Chris Gorgolewski 2012

classdef NegativeGamma < Gamma

    methods
        function obj = NegativeGamma(shape, scale, mu)
            if nargin == 0
                obj.shape = 5;
                obj.scale = 10;
                obj.mu = 0;
            else
                obj.shape = shape;
                obj.scale = scale;
                obj.mu = mu;
            end
        end
        function fit_weighted(obj, data, weights)
            [obj.shape, obj.scale] = gam_param(-(data - obj.mu), weights);
        end
        function y = pdf(obj, data)
            y = gam_dens(obj.shape, obj.scale, -(data - obj.mu));
        end
    end
    
end

