function [ W ] = build_mdas_W(ana_par, fKappa, FKappa)
%BUILD_MDAS_W Summary of this function goes here
%   Detailed explanation goes here

    K0 = ana_par.K_max;
    Keps = ana_par.kappa_eps;
    
    W_raw = @(y) fKappa(K0 - y)./(1 - FKappa(K0 - y));
    
    ulim = ana_par.K_max - ana_par.K_min;
    llim = 0;
    KK = linspace(llim, ulim, 80*ana_par.n_K);
    
    WW = W_raw(KK);
    
    KK1 = [llim - 10*Keps, llim - Keps, KK, ulim + Keps, ulim + 10*Keps];
    WW1 = [0, 0, WW, 0, 0];
    
    W = @(y) interp1(KK1, WW1, y, 'linear', 'extrap');
end

