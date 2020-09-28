function [kappa, kappa_inv, kappa_prime, kappa_prime_inv] = build_kappas(par, e_inv_asc)
%
% define the kappa functions, necessary for inverting the relationship
% between K^+_a and sigma_a
%

    M_qui = @(s) (s < par.s_c) ;
    kappa_exact = @(s) s./e_inv_asc(s) .* M_qui(s);
    
    %S = par.S;
    S = linspace(0, par.s_c, 20*par.n_s);
    S = S(2:(end-1));
    K = kappa_exact(S);
    
    S_aug = [S, par.s_c + 1e-3, 10*par.s_c];
    K_aug = [K, 1, 1];
    
    % use makima interpolation, because we want a continuous first
    % derivative later
    kappa = @(s) interp1(S_aug, K_aug, s, 'makima', 'extrap');
    
    
    
    
    K_eps = kappa(0);
    par.K_max = K_eps;
    par.K_min = kappa(par.s_c - par.kappa_eps);
    par.update_derived_parameters();
    
    
    
    
    
    S_aug = [S, par.s_c, par.s_c];
    K_aug = [K, 1 - 1e-3, 0];
    
    kappa_inv = @(k) interp1(K_aug, S_aug, k, 'makima', 'extrap');
    
    
    
    DK = diff(K);
    DS = diff(S);
    muS = (S(2:end) + S(1:(end-1)))/2;
    muK = (K(2:end) + K(1:(end-1)))/2;
    
    KpS = DK./DS;
    
    KpS_aug =  [-1, -1, KpS, -1e6, -1e7];
    muS_aug = [-1, 0, muS, 1, 2];
    muK_aug = [10, par.K_max, muK, 1, 0];
    
    kappa_prime = @(s) interp1(muS_aug, KpS_aug, s, 'makima', 'extrap');
    kappa_prime_inv = @(k) interp1(muK_aug, KpS_aug, k, 'makima', 'extrap');
    
    
    % alternate methods were bad
    
    %switch par.kappa_inversion_rule
    %end
    
end