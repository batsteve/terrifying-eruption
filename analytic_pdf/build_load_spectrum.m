function [fS, fS_qui, fS_ex, sampler] = build_load_spectrum(par)
%
% Define the load spectrum
%
% probability density of a spike of magnitude s (with tolerance ds) per 
% unit time
%

    sampler = @(D) ones(D);

    M_qui = @(s) (s < par.s_c);

    
    switch par.load_spectrum
            
        case 'comparison-independent-spike-model'
            
            p = 1/20000;
            mu = 0.03;
            g = 0.03;
            r = 0.42;
            
            fS_raw = @(s) (1-p)*normpdf(s, mu, g)  + ...
                (1-p)*normpdf(-s, mu, g) + ...
                + p*raylpdf(s, r);
            
            sampler = @(D) sample_func(D, mu, g, p, r);
            
            

        otherwise
            warning('%s not recognized!\n');
    end
    

    fS_norm = integrate_spectrum(par, fS_raw, par.s_over_factor.*par.s_c);
    
    fS = @(s) fS_raw(s)/fS_norm;
    
    
    %
    % construct conditional pdfs, given that there is no extreme spike
    %
    
    C_qui = integrate_spectrum(par, @(s) fS(s).*M_qui(s), par.s_c);

    fS_qui = @(s) fS(s)./C_qui.*M_qui(s);
    fS_ex = @(s) fS(s)./(1-C_qui).*(1 - M_qui(s));
    
end

function [ yy ] = sample_func(D, a1, a2, a3, a4)
    Q = a1 + a2*randn(D);
    Q = max(abs(Q), 1e-4);
    EX = binornd(1, a3, D);
    X =  raylrnd(a4, D);

    yy = (1 - EX).*Q + EX.*X;
end

