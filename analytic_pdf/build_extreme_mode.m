function [f_ex] = build_extreme_mode(par, fS, M_qui)
%BUILD_EXTREME_MODE Summary of this function goes here
%   Detailed explanation goes here
    
    %f_ex = @(t) C_qui.^(t/dt-1).*C_qui;
    
    P_qui_1 = integrate_spectrum(par, @(s) fS(s).*M_qui(s), par.s_c);
    %P_qui_2 = integral(fS, par.s_c, inf);
    
    lambda = (1 - P_qui_1);
    f_ex = @(t) lambda*exp(-lambda * t);

end

