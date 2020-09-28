function [DS] = build_independent_spike_signal(so_par, seed, sampler)
%BUILD_INDEPENDENT_SPIKE_SIGNAL Summary of this function goes here
%   Detailed explanation goes here

    if (~isnan(seed))
        rng(seed);
    end
    
    n = ceil(so_par.input_signal_length);
    
    yy = sampler([n, 1]);
    
    DS = struct;
    
    DS.sspp = yy';
    DS.ssnp = zeros(size(yy));
    
end

