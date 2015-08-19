% EM algorithm, bic and plot
% ------------------------------------------
% This file is part of the adaptative_threshold SPM plug-in
% Copyright (C) Cyril Pernet and Chris Gorgolewski 2012

classdef ExpectationMaximization < handle
    properties
        components
        mix
    end
    methods
        
        function obj = ExpectationMaximization(components)
            % update the model component
            obj.components = components;
            obj.mix = ones(1,length(obj.components)) * 1 / length(obj.components);
        end
        
        function posteriors = posteriors(obj, data)
            % update the model posterior prob. of the data given the model
            posteriors = zeros(length(data), length(obj.components));
            for i=1:length(obj.components)
                posteriors(:, i) = obj.components{i}.pdf(data);
            end
            posteriors(isnan(posteriors)) = 0;
            posteriors = posteriors .* repmat(obj.mix, length(posteriors), 1);
        end
        
        function ll = loglikelihood(obj, data)
            % likelihood of the posterior
            likelihood = obj.posteriors(data);
            sum_likelihood = sum(likelihood,2);
        	ll = sum(log(sum_likelihood));
        end

        function bic = BIC(obj, data)
            % Bayesian information criterion
            k = 0;
            for i=1:length(obj.components)
                k = k + obj.components{i}.free_parameters;
            end
            bic = -2*obj.loglikelihood(data) + k * log(length(data));
        end

        function fit(obj, data)
            % fit the model to the data
            maxiter=100;
            min_improvement=1.e-4;
            prev_loglike = 0;

            for i=1:maxiter
 
                resp = obj.E(data);
                obj.M(data, resp);

                cur_loglike = obj.loglikelihood(data);
                improv = cur_loglike - prev_loglike;
                if i > 1 && improv <= min_improvement
                    fprintf('break %f\n', improv)
                    break
                end
                prev_loglike = cur_loglike;
                fprintf('log likelihood = %f, improv = %f, iter = %d\n',obj.loglikelihood(data), improv, i)
            end
        end
        
        function plot(obj, data)
            % plot the data (histrogram) and the model
            hold all; histnorm(data,50); colormap([.5 .5 1]);
            x = linspace(min(data), max(data), 1000);
            pdf_sum = zeros(1, length(x));
            for i=1:length(obj.components)
                 if i == 1; c = 'r'; elseif i ==2; c ='g'; else  c = 'y'; end
                plot(x, obj.components{i}.pdf(x).*obj.mix(i),c,'Linewidth',3);
                pdf_sum = pdf_sum + obj.components{i}.pdf(x).*obj.mix(i);
            end
            plot(x, pdf_sum,'--k','LineWidth',2); 
            
            if size(obj.components,2) == 1
                mytitle = sprintf('Gaussian only \n BIC=%g', obj.BIC(data));
            elseif size(obj.components,2) == 2
                mytitle = sprintf('Gaussian & Gamma \n BIC=%g', obj.BIC(data));
            elseif size(obj.components,2) == 3
                mytitle = sprintf('Gaussian & +/-Gamma \n BIC=%g', obj.BIC(data));
            end
            title(mytitle,'Fontsize',16);
            grid on; hold off
        end
    end
    
    methods(Access=private)
        function resp = E(obj, data)
            resp = obj.posteriors(data);
            resp_sum = sum(resp,2);
            resp = resp ./ repmat(resp_sum, 1, length(obj.components));
        end
        
        function M(obj, data, resp)
            for i=1:length(obj.components)
                if (strcmp(class(obj.components{i}), 'NegativeGamma') || strcmp(class(obj.components{i}), 'Gamma'))
                    for j=1:length(obj.components)
                        if strcmp(class(obj.components{j}), 'Gaussian')
                            obj.components{i}.mu = obj.components{j}.mu;
                            break
                        end
                    end
                end
                obj.components{i}.fit_weighted(data, resp(:, i))
                obj.mix(i) = sum(resp(:, i)) / length(resp(:, i));
            end
        end
    end
end



    