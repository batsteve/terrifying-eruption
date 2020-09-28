function [ f_nfail ] = build_quiescent_mode_hvqs(ana_par, Feta, V, Vexp, Q_prime, S)
%BUILD_QUIESCENT_MODE_VQS Summary of this function goes here
%   Detailed explanation goes here

    dk = abs(ana_par.delta_K);
    K0 = ana_par.K_max;

    switch ana_par.Veta_exp_rule
            %integrand = @(xi, zeta, n) Feta(K0 - xi - n*dk).*V(zeta + xi, n*dk + xi).*Q_prime(xi).*S(zeta);
        case 'taylor-expansion'
            integrand = @(xi, zeta, n) (1 - Feta(K0 - xi - n*dk)).*(1 + V(zeta + xi, n*dk + xi)).*Q_prime(xi).*S(zeta);
    
        case 'exponential'
            integrand = @(xi, zeta, n) (1 - Feta(K0 - xi - n*dk)).*Vexp(zeta + xi, n*dk + xi).*Q_prime(xi).*S(zeta);
            
        otherwise
            warning('%s not recognized!\n', ana_par.Veta_exp_rule)
    end
    
    
    switch ana_par.hvqs_lims
        case 'zero-k-max'
            xi_ulim = ana_par.K_max;
            xi_llim = 0;
            zeta_ulim = ana_par.K_max;
            zeta_llim = 0;
        case 'infinity'
            xi_ulim = inf;
            xi_llim = -inf;
            zeta_ulim = inf;
            zeta_llim = -inf;
        case 'zero-phi'
            xi_ulim = ana_par.K_max - ana_par.K_min;
            xi_llim = 0;
            zeta_ulim = ana_par.K_max - ana_par.K_min;
            zeta_llim = 0;
        otherwise
    end
    
    f_n_exact = @(n) 1/dk*integral2( @(xi, zeta) integrand(xi, zeta, n), ...
        xi_llim, xi_ulim, zeta_llim, zeta_ulim, ...
        'AbsTol', ana_par.f_qui_abs_tol, 'RelTol', ana_par.f_qui_rel_tol);
    
    NN = linspace(0, ana_par.T_max, ana_par.n_t);
    PP = zeros(size(NN));
    for k = 1:length(NN)
        tic;
        fprintf('    Starting t step %d (out of %d).', k, ana_par.n_t);
        PP(k) = f_n_exact(NN(k));
        fprintf('  (%0.2fs).\n', toc);    
    end
    
    dt = NN(2) - NN(1);
    C = sum(PP(:))*dt;
    
    PP_norm = (1/C)*PP;
    
    p = 0; %extrap value
    f_nfail = @(t) interp1(NN, PP_norm, t, 'makima', p);
    
    %X = linspace(0, ana_par.K_max, ana_par.n_K);
    %Z = linspace(0, ana_par.K_max, ana_par.n_K);
    %
    %[ f_nfail ] = build_interp2(f_n_exact, X, Z, 0);
end

