function [fEta, FEta] = build_eta_pdf(ana_par, fS, eta, eta_inv, eta_prime, eta_prime_inv)
%BUILD_KAPPA_PDF Summary of this function goes here
%   Detailed explanation goes here

    
    
    switch ana_par.kappa_eta_pdf_formula
        case 'prime-inv-combined'
            pdf_expression = @(e) fS(eta_inv(e))./abs(eta_prime_inv(e));
            
        case 'prime-inv-separate'
            pdf_expression = @(e) fS(eta_inv(e))./abs(eta_prime(eta_inv(e)));
        otherwise
    end
    
    
    E = linspace(0, ana_par.K_min, 20*ana_par.n_K);
    de = E(2) - E(1);
    
    fE = pdf_expression(E);
                

    FE = cumsum(fE)*de;
    
    C = FE(end);
    FE_adj = (1/C) * FE;
    fE_adj = (1/C) * fE;
    
    p = 0; % extrap value for a distribution vanishes
    fEta = @(k) interp1(E, fE_adj, k, 'makima', p);
    
    Keps = ana_par.kappa_eps;
    FE_aug = [0, 0, FE_adj(2:end), 1, 1];
    E2 = [-1, 0, E(2:end), ana_par.K_min + Keps, 5*ana_par.K_min];
    
    FEta = @(e) interp1(E2, FE_aug, e, 'linear', 'extrap');
    
end

