% This code was ported from the nipy project (http://nipy.org) and was
% originally contrributed by Bertrand Thirion - returns gamma parameters

function [shape, scale] = gam_param(data, weights)
    eps = 1.e-5;
    i = data > 0;
    szi = sum(weights(i));
    if szi > 0
        shape = compute_c(data(i), weights(i), eps);
        scale = dot(data(i),weights(i)) / (szi * shape);
    else
        shape = 1;
        scale = 1;
    end
end

function scale = compute_c(data, weights, eps)
    y = dot(weights,log(data)) / sum(weights) - log(dot(weights,data) / sum(weights));
    eps = 1.e-7;
    if y > - eps
        scale = 10;
    else
        scale = psi_solve(y, 0.00001);
    end
end
    
function r = psi_solve(y, eps)
    if y > 0
        throw(MException('y>0, the problem cannot be solved'))
    end
    u = 1;
    if y > psi(u) - log(u)
        while psi(u) - log(u) < y
            u = u * 2;
        end
        u = u / 2;
    else
        while psi(u) - log(u) > y
            u = u / 2;
        end
    end
    r = dichopsi_log(u, 2 * u, y, eps);
end
    
function t = dichopsi_log(u, v, y, eps)
    if u > v
        tmp = u;
        u = v;
        v = tmp;
    end
    t = (u + v) / 2;
    if abs(u - v) >= eps
        if psi(t) - log(t) > y
            t = dichopsi_log(u, t, y, eps);
        else
            t = dichopsi_log(t, v, y, eps);
        end
    end
end