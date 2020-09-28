function [ q_prod, q_prod_prime, q_integrand ] = build_mdas_q_prod(ana_par, FKappa)
%BUILD_MDAS_Q_PROD Summary of this function goes here
%   Detailed explanation goes here

    dk = ana_par.delta_K;
    K0 = ana_par.K_max;
    K1 = ana_par.K_min;
    
    ulim = K0 - K1; % this is the most that phi can be
                    % should be equal to phi_max * delta_K
    
    integrand_exact = @(z) real(log(1 - FKappa(K0 - z)));
    
    ZZ = linspace(0, ulim, 80*ana_par.n_K);
    II = zeros(size(ZZ));
    for k = 1:length(ZZ)
        II(k) = integrand_exact(ZZ(k));
    end
    
    ZZ1 = [ZZ, ulim + dk, ulim + 10*dk];
    II1 = [II, 0, 0];
    
    q_integrand = @(y) interp1(ZZ1, II1, y, 'makima', 'extrap');
    
    
    
    q_prod_raw = @(y) exp(1/dk*integral(q_integrand, y, ulim));
    
    KK = linspace(0, ulim, 20*ana_par.n_K);
    QQ = zeros(size(KK));
    for k = 1:length(KK)
        QQ(k) = q_prod_raw(KK(k));
    end
    
    KK1 = [-1, KK, ulim + 1];
    QQ1 = [0, QQ, 1];
    
    q_prod = @(y) interp1(KK1, QQ1, y, 'makima', 'extrap');
    
    
    
    switch ana_par.q_prod_prime_algorithm
        case 'finite-difference'
            Keps = 2e-2;
            QQp = zeros(size(KK));
            for k = 2:(length(KK)-1)
                QQp(k) = 1/(2*Keps)*(q_prod_raw(KK(k) + Keps) - q_prod_raw(KK(k) - Keps));
            end
            QQp(1) = 1/(Keps)*(q_prod_raw(KK(1) + Keps) - q_prod_raw(KK(1)));
            QQp(end) = 1/(Keps)*(q_prod_raw(KK(end)) - q_prod_raw(KK(end) - Keps));

            p = 0; % extrap value
            q_prod_prime = @(y) interp1(KK, QQp, y, 'linear', p);
    
        case 'analytic'
            derviative_exact = @(y) -1/dk*q_prod(y).*q_integrand(y);
            
            QQp = zeros(size(KK));
            for k = 2:(length(KK)-1)
                QQp(k) = derviative_exact(KK(k));
            end
            
            % catch weird things?
            QQp(QQp < 0) = 0;
            
            p = 0; % extrap value
            q_prod_prime = @(y) interp1(KK, QQp, y, 'linear', p);
            
        otherwise
            warning('%s not recognized!\n', ana_par.q_prod_prime_algorithm);
    end
    
end

