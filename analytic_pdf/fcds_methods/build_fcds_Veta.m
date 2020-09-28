function [ Veta, integrand, Vexp ] = build_fcds_Veta(ana_par, Feta)
%BUILD_FCDS_VETA Summary of this function goes here
%   Detailed explanation goes here

    dk = ana_par.delta_K;
    K0 = ana_par.K_max;
    Keps = ana_par.kappa_eps;
    
    integrand_exact = @(s) log(Feta(K0 - s));
    
    
    
    switch ana_par.Veta_interp_region
        case 'zero-k-max'
            a1 = 0;
            a2 = ana_par.K_max;
        case 'min-k-max'
            a1 = ana_par.K_min;
            a2 = ana_par.K_max;
        case 'phi-k-max'
            a1 = ana_par.K_max - ana_par.K_min;
            a2 = ana_par.K_max;
        case 'fat-phi-k-max'
            delta = (ana_par.K_max - ana_par.K_min) / 20;
            a1 = ana_par.K_max - ana_par.K_min - delta;
            a2 = ana_par.K_max + delta;
            
        otherwise
    end
    
    
    switch ana_par.Veta_integrand_rule
        case 'exact'
            integrand = @(s) integrand_exact(s);
            
        case 'interp'
            
            S = linspace(a1, a2, 80*ana_par.n_K);
            I = integrand_exact(S);
            
            S1 = [a1 - 10*Keps, a1 - Keps, S, a2 + Keps, a2 + 10*Keps];
            %I1 = [-inf, -100, I, 0, 0];
            I1 = [0, 0, I, -100, -100];
            
            I1 = I1(isfinite(I1));
            S1 = S1(isfinite(I1));
            
            %p = -1;  % extrap value
            integrand = @(s) interp1(S1, I1, s, 'linear', 'extrap');
            
            
        otherwise
            warning('%s not recognized!\n', ana_par.Veta_integrand_rule)
    end
            
            
            
            
    %Veta_exact_1 = @(xi, n, na) integral(integrand, na*dk+xi, n*dk+xi);
    Veta_exact_2 = @(llim, ulim) 1/dk*integral(integrand, llim, ulim);
    
    Vexp_exact = @(llim, ulim) exp(1/dk*integral(integrand, llim, ulim));
    
    
    
    switch ana_par.Veta_interp_region
        case 'zero-k-max'
            XX = linspace(0, ana_par.K_max, ana_par.n_K);
        case 'min-k-max'
            XX = linspace(ana_par.K_min, ana_par.K_max, ana_par.n_K);
        case 'phi-k-max'
            XX = linspace(ana_par.K_max - ana_par.K_min, ana_par.K_max, ana_par.n_K);
            
        otherwise
    end
            
    

    switch ana_par.Veta_algorithm
        case 'exact'
            Veta = @(x1, x2) Veta_exact_2(x1, x2);
            Vexp = @(x1, x2) Vexp_exact(x1, x2);
            
        case 'interp'
            [XX1, XX2] = meshgrid(XX, XX);
            VV = zeros(length(XX), length(XX));
            
            % figure we only need the upper triangle of this array?
            % function should be skew-symmetric
            
            for k1 = 1:length(XX)       
                for k2 = k1:length(XX)
                    VV(k1, k2) = Veta_exact_2(XX1(k1, k2), XX2(k1, k2));
                    %VV(k2, k1) = -VV(k1, k2);  % should backwards integral change sign?
                    VV(k2, k1) = 0;  % should backwards integral change sign?
                end
            end
            
            %Veta = @(x1, x2) build_interp2(f, X, Y, p)
            p = 0;  % what should this actually be?
            Veta = @(x1, x2) interp2(XX1,XX2,VV,x1,x2, 'makima');
            % also compute exp(V)
            % because it might be important
            
            
            for k1 = 1:length(XX)
                for k2 = k1:length(XX)
                    VV(k1, k2) = Vexp_exact(XX1(k1, k2), XX2(k1, k2));
                    VV(k2, k1) = 0;   % Solves some problems, maybe?
                end
            end

            %Veta = @(x1, x2) build_interp2(f, X, Y, p)
            p = 0;  % what should this actually be?
            Vexp = @(x1, x2) interp2(XX1,XX2,VV,x1,x2, 'makima');
            
        case 'sum-of-antiderivatives'  
            bb = (1/2)*(max(S) - min(S));
            
            V_antiderivative_exact = @(x) Veta_exact_2(bb, x);
            
            VV = zeros(size(S));
            for k = 1:length(S)
                VV(k) = V_antiderivative_exact(S(k));
            end
            
            V_antid = @(x) interp1(S, VV, x, 'makima', 'extrap');
            
            Veta = @(x1, x2) min(V_antid(x2) - V_antid(x1), zeros(size(x1)));
            Vexp = @(x1, x2) min(exp(V_antid(x2) - V_antid(x1)), ones(size(x1)));
            %Vexp = @(x1, x2) exp(V_antid(x2) - V_antid(x1));
    end
    
  
    % not quite sure what to do about integration bounds that are bigger
    % (or smaller!) than the sensible support of the integrand.  This idea
    % doesn't seem to be the right choice
    %
    % But after a very large amount of senseless bugchasing, I'm pretty
    % sure that we do need some kind of clip to avoid the case where
    % (x1 > a2).  Why this and only this is the problem, I don't know.  Why
    % it took be two days to find this, I don't know.
    switch ana_par.Veta_exp_clip
        case true
            % Bad ideas!  No cookie!
            %Veta = @(x1, x2) Veta(x1, x2).*(x1 >= a1).*(x2 <= a2);
            %Vexp = @(x1, x2) Vexp(x1, x2).*(x1 >= a1).*(x2 <= a2) + 1.*~((x1 >= a1).*(x2 <= a2));
            %Vexp = @(x1, x2) 1 + (1 - Vexp(x1, x2)).*(x1 >= a1).*(x2 <= a2).*(x2 >= x1);
            %Vexp = @(x1, x2) Vexp(x1, x2).*(x2 >= x1);
            
            
            
            Vexp = @(x1, x2) Vexp(x1, x2).*(x2 <= a2);
        otherwise
            % do nothing
    end
    
    
end

