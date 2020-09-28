function [ Szeta ] = build_fcds_Szeta( ana_par, Q, W)
%BUILD_FCDS_SZETA Summary of this function goes here
%   Detailed explanation goes here

    switch ana_par.szeta_lims
        case 'infinity'
            ulim = inf;
            llim = -inf;
            
        case 'zero-k-max'
            ulim = ana_par.K_max;
            llim = 0;
            
        otherwise
            warning('%s not recognized!\n', ana_par.szeta_lims);
    end

    integrand = @(zeta, y) Q(y).*W(zeta+y);
    S_exact = @(zeta) integral(@(y) integrand(zeta, y), llim, ulim);
    
    ZZ = linspace(0, ana_par.K_max, 20*ana_par.n_K);
    %SS = S_exact(ZZ);
    SS = zeros(size(ZZ));
    for k = 1:length(ZZ)
        SS(k) = S_exact(ZZ(k));
    end
    
    
    %p = 0; % extrap value
    Szeta = @(zeta) interp1(ZZ, SS, zeta, 'linear', 'extrap');
    
end

