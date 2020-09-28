function [delta_K] = estimate_delta_k(par, fS)
%ESTIMATE_DELTA_K Summary of this function goes here
%   Detailed explanation goes here

    s_bar = integral(@(s) s.*fS(s), 0, inf, 'ArrayValued', 1);
    
    delta_K = s_bar/par.d_a;

end

