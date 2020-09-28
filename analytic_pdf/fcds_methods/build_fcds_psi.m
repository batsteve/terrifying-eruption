function [ psi, nDS ] = build_fcds_psi(ana_par, e_inv_des)
%BUILD_FCDS_PSI Summary of this function goes here
%   Detailed explanation goes here

    delta_K = ana_par.delta_K;
    
    %
    % Define the psi function, which calculates the anticipation
    %

    psi = @(s, delta_n, k) s - (k - delta_n'*delta_K).*e_inv_des(s);

    %
    % partial inversion
    %

    nDS = @(s, k) (k - s./e_inv_des(s))./(delta_K);
end

