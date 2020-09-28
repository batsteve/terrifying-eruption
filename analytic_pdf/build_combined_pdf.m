function [ f_tot ] = build_combined_pdf( ana_par, f_qui, f_ex )
%BUILD_COMBINED_PDF Summary of this function goes here
%   Detailed explanation goes here

    TT = linspace(0, ana_par.T_max, ana_par.n_t);
    dt = TT(2) - TT(1);

    F_qui = f_qui(TT)';

    F_qui = max(F_qui, 0);
    C_qui_2 = sum(F_qui*dt);
    F_qui = F_qui/C_qui_2;

    
    G_qui = cumsum(F_qui)*dt;

    F_ex = f_ex(TT)';
    F_ex = max(F_ex, 0);
    G_ex = cumsum(F_ex)*dt;
    
    F_full = F_qui.*(1-G_ex) + F_ex.*(1 - G_qui);
    
    f_tot = @(y) interp1(TT,F_full,y, 'linear', 'extrap'); 
    
    
end

