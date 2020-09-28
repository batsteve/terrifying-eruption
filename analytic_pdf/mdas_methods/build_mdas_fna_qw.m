function [ fna ] = build_mdas_fna_qw(ana_par, Q, W )
%BUILD_MDAS_PHI_QW Summary of this function goes here
%   Detailed explanation goes here

    dk = ana_par.delta_K;

    %f_na_raw = @(a) integral(@(x) Q(x * dk) .* W((a + x).*dk).*dk, -inf, inf);
    
    ulim = ana_par.K_max - ana_par.K_min;
    f_na_raw = @(a) 1/dk*integral(@(xi) Q(xi) .* W(a.*dk + xi).*dk, -0, ulim);
    
    NN = linspace(0, ana_par.T_max, 10*ana_par.n_t);
    %PP = f_phi_raw(NN);
    PP = zeros(size(NN));
    for k = 1:length(NN)
        PP(k) = f_na_raw(NN(k));
    end
    
    fna = @(y) interp1(NN, PP, y, 'makima', 'extrap');
    
    
end

