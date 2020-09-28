function [ fS_bin ] = max_bin_spectrum(par, fS_per_spike)
%MAX_BIN_SPECTRUM Summary of this function goes here
%   Detailed explanation goes here

    %fS_per_spike = @(s) 2*(1-p)*normpdf(s, 0, g) + p*raylpdf(s, r);
    
    
    gS_per_spike = @(s) integral(fS_per_spike, 0, s, 'ArrayValued', 1);
    
    n_bin = par.T_max/par.n_t;
    fS_exact = @(s) n_bin*(gS_per_spike(s).^n_bin)*fS_per_spike(s);
    
    Sf = par.S;
    ff = zeros(length(Sf), 1);
    for k = 1:length(Sf)
        ff(k) = fS_exact(Sf(k));
    end
    fS_bin_raw = @(s) interp1(Sf,ff,s, 'linear', 'extrap');
    
    [C] = integrate_spectrum(par, fS_bin_raw, par.s_over_factor.*par.s_c);
    
    fS_bin = @(s) fS_bin_raw(s)/C;
    
end

