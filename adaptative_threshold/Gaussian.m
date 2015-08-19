% create Gaussian model
% ------------------------------------------
% This file is part of the adaptative_threshold SPM plug-in
% Copyright (C) Cyril Pernet and Chris Gorgolewski 2012

classdef Gaussian < Distribution
    properties
        mu = 0
        sigma = 1
        sigma_sqr = 1
    end
    properties(Constant)
        free_parameters = 2
    end
    methods
        function obj = Gaussian(mu,sigma)
            obj.mu = mu;
            obj.sigma = sigma;
            obj.sigma_sqr = sigma * sigma;
        end
        function y = pdf(obj, data)
            p1 = -.5 * (power(data - obj.mu,2) / obj.sigma_sqr);
            p2 = (obj.sigma * sqrt(2*pi));
            y = exp(p1) ./ p2; 
        end
        function fit(obj, data)
            obj.mu = sum(data)/length(data);
            obj.sigma = sqrt(sum((data - obj.mu).^2)/length(data));
        end
        function fit_weighted(obj, data, weights)
            weights_sum = sum(weights);
            obj.mu = sum(data .* weights) / weights_sum;
            obj.sigma_sqr = sum(power((data - obj.mu),2) .* weights) / weights_sum;
            obj.sigma = sqrt(obj.sigma_sqr);
        end
    end
end


