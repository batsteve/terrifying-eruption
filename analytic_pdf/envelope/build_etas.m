function [eta, eta_inv, eta_prime, eta_prime_inv] = build_etas(ana_par, e_inv_des)
%BUILD_FCDS_ETA Summary of this function goes here
%   Detailed explanation goes here

    M_qui = @(s) (s < ana_par.s_c) ;
    eta_exact = @(s) s./e_inv_des(s) .* M_qui(s);
    
    S = linspace(0, ana_par.s_c, 20*ana_par.n_s);
    S = S(2:(end-1));
    E = eta_exact(S);
    
    S_aug = [-1, -1e-3, S, ana_par.s_c, 2*ana_par.s_c];
    E_aug = [0, 0, E, 1, 1];
    
    % use makima interpolation, because we want a continuous first
    % derivative later
    eta = @(s) interp1(S_aug, E_aug, s, 'makima', 'extrap');

    S_aug = [0, 0, S, ana_par.s_c, ana_par.s_c];
    E_aug = [-10, -1, E, 1, 10];
    
    eta_inv = @(e) interp1(E_aug, S_aug, e, 'makima', 'extrap');

    
    DE = diff(E);
    DS = diff(S);
    muS = (S(2:end) + S(1:(end-1)))/2;
    muE = (E(2:end) + E(1:(end-1)))/2;
    
    EpS = DE./DS;
    
    EpS_aug = [EpS];
    muS_aug = [muS];
    muE_aug = [muE];
    
    eta_prime = @(s) interp1(muS_aug, EpS_aug, s, 'linear', 'extrap');
    eta_prime_inv = @(e) interp1(muE_aug, EpS_aug, e, 'linear', 'extrap');
    
    % other methods were bad
    
    %switch ana_par.eta_prime_fd_algorithm
    %end
    
    
    
    
    
%     switch ana_par.eta_prime_fd_algorithm
%         case 'method-1'
%     
% 
%             Sp = linspace(5*ana_par.ds, ana_par.s_c - ana_par.ds, 10*ana_par.n_s);
% 
%             EpS = zeros(length(Sp), 1);
%             for k = 1:length(Sp)
%                 EpS(k) = eta_prime_exact(Sp(k), eta, ana_par);
%             end
%             eta_prime = @(s) interp1(Sp, EpS, s, 'linear', 'extrap');
%             
%         case 'method-2'
%             MuS = 1/2*(S(1:(end-1)) + S(2:end));
%             DS = diff(S);
%             DES = diff(ES);
%             
%             eta_prime = @(s) interp1(MuS, DES./DS, s, 'linear', 'extrap');
%             
%         otherwise
%             warning('%s not recognzied!\n', ana_par.eta_prime_fd_algorithm);
%     end
%     
%     
%     
%     
%     
%     
%     
%     Epi = eta_prime(eta_inv(E2));
%     
%     switch length(E2)
%         case 1
%             eta_prime_inv = @(e) eta_prime(0.1);
%         otherwise
%             eta_prime_inv = @(e) interp1(E2, Epi, e, 'linear', 'extrap');
%     end
%     
end

% function [ kp ] = eta_prime_exact(s, eta, ana_par)
%     ds = min([s/2; ana_par.s_c - s/2; ana_par.ds]);
%     kp = (eta(s + ds) - eta(s - ds))/(2*ds);
% end