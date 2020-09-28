function [fKappa, FKappa] = build_kappa_pdf(ana_par, fS, kappa, kappa_inv, kappa_prime, kappa_prime_inv)
%BUILD_KAPPA_PDF Summary of this function goes here
%   Detailed explanation goes here


    switch ana_par.kappa_eta_pdf_formula
        case 'prime-inv-combined'
            pdf_expression = @(k) fS(kappa_inv(k))./abs(kappa_prime_inv(k));
            
        case 'prime-inv-separate'
            pdf_expression = @(k) fS(kappa_inv(k))./abs(kappa_prime(kappa_inv(e)));
        otherwise
    end
    
    K = linspace(ana_par.K_min, ana_par.K_max, 20*ana_par.n_K);
    dk = K(2) - K(1);
    fK = pdf_expression(K);
                
    p = 0; % extrap value for a distribution vanishes
    fKappa = @(k) interp1(K, fK, k, 'makima', p);
    
    FK = cumsum(fK)*dk;
    
    Keps = ana_par.kappa_eps;
    FK = [0, 0, FK, 1, 1];
    K2 = [0, ana_par.K_min - Keps, K, ana_par.K_max + Keps, 5*ana_par.K_max];
    
    FKappa = @(k) interp1(K2, FK, k, 'linear', 'extrap');
    
end

