function [phi, phi_inv, phi_prime, phi_inv_theta] = ...
    build_mdas_phi(par, kappa, kappa_inv, kappa_prime)
% Define the phi function for leading spike damage
%

    delta_K = par.delta_K;


    K_eps = kappa(0);

    %phi = @(s, t) K_eps/delta_K - s ./ (e_inv_asc(s) * delta_K) - t;
    phi = @(s, t) (K_eps - kappa(s)) / delta_K - t;

    phi_max = phi(par.s_c - par.sigma_eps, 0);        % useful for graphing
    

    
    phi_inv = @(p, t) kappa_inv(max(K_eps - (p + t)*delta_K, 0));
    phi_inv_theta = @(theta) kappa_inv(K_eps - theta*delta_K );
    
    
    
    
    phi_prime = @(s, t) - kappa_prime(s)/ delta_K;

    
    
    par.phi_max = phi_max;
    par.update_derived_parameters();
end

