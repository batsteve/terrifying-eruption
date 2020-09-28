function [ results ] = build_full_pdf(par, fS, e1, e2)
%BUILD_FULL_PDF Computes the failure time pdf for the Serebrinsky-Ortiz
%model of material fatigue failure
%   par -- parameter object
%   fS -- function representing the distribution of load peaks
%   either:
%       1) e1 -- functional representation of the coherent envelope in the
%       form e(delta)
%       2) e1, e2 -- functional representions of the inverted ascending and
%       descending legs of the coherent envelope, in the form e(sigma)

    % quiescent spike mask (kills extreme spikes)
    M_qui = @(s) (s < par.s_c);
    
    C_qui = integrate_spectrum(par, @(s) fS(s).*M_qui(s), par.s_c);
    fS_qui = @(s) fS(s)./C_qui.*M_qui(s);

    if (par.plot_verbosity >= 1)
        Plotter.draw_spectrum_pdf(fS, par.S);
        %Plotter.draw_spectrum_cdf(fS, par.S);
    end
    
    
    switch nargin
        case 3
            fprintf('Inverting coherent envelope.\n');
            e_forward = e1;
            [ e_inv_asc, e_inv_des ] = invert_envelope( par, e_forward );
            
        case 4
            e_inv_asc = e1;
            e_inv_des = e2;

        otherwise
            warning('Wrong number of input arguments!  (Expected 3 or 4, got %d).\n', nargin);
    end
    
    if (par.plot_verbosity >= 2)
        Plotter.draw_coherent_envelope(e_inv_asc, e_inv_des, par.s_c, par.d_c);
    end
    
    
    
    

    fprintf('Building kappas and etas.\n');
    [kappa, kappa_inv, kappa_prime, kappa_prime_inv] = build_kappas(par, e_inv_asc);
    [fKappa, FKappa ] = build_kappa_pdf(par, fS_qui, kappa, kappa_inv, kappa_prime, kappa_prime_inv);
    
    [eta, eta_inv, eta_prime, eta_prime_inv] = build_etas(par, e_inv_des);
    [fEta, FEta ] = build_eta_pdf(par, fS_qui, eta, eta_inv, eta_prime, eta_prime_inv);

    
    if (par.plot_verbosity >= 1)
        Plotter.draw_kappa_pdf(fKappa, FKappa, par.K);
        Plotter.draw_eta_pdf(fEta, FEta);
    end
    
    if (par.plot_verbosity >= 2)
        %Plotter.draw_kappa(kappa, par.S);
        Plotter.draw_kappa_family(kappa, kappa_inv, kappa_prime, kappa_prime_inv, par.S, par.K);
        Plotter.draw_eta_family(eta, eta_inv, eta_prime, eta_prime_inv, par.S);
        Plotter.draw_envelope_kappa_eta_transformation(kappa, eta, par.s_c);
        %Plotter.draw_envelope_kappa_eta_comparison(e_inv_asc,e_inv_des, kappa, eta, par.s_c);
    end


    fprintf('Computing MDAS quantities.\n')

    tic;
    fprintf('  Computing phi function.')
    [phi, phi_inv, phi_prime, phi_inv_theta] = build_mdas_phi(par, kappa, kappa_inv, kappa_prime);
    fprintf('  (%0.2fs).\n', toc);
    
    tic;
    fprintf('  Computing Q_prod function.')
    [ q_prod, q_prod_prime, q_integrand ] = build_mdas_q_prod(par, FKappa);
    fprintf('  (%0.2fs).\n', toc);
        
    tic;
    fprintf('  Computing W function.');
    [ W ] = build_mdas_W(par, fKappa, FKappa);
    fprintf('  (%0.2fs).\n', toc);
    
    tic;
    fprintf('  Computing fna distribtion.');
    [ fna ] = build_mdas_fna_qw(par, q_prod, W );
    fprintf('  (%0.2fs).\n', toc);
    
    
    
    
    if (par.plot_verbosity >= 1)
        Plotter.draw_phi_truncated(phi, par.S, par.T, M_qui);
        Plotter.draw_q_prod(q_prod, par.K_max- 1);
        Plotter.draw_fna(fna, par.T);
    end
    if (par.plot_verbosity >= 2)

        Plotter.draw_q_integrand(q_integrand, par.K_max - 1);
        Plotter.draw_q_prod_prime(q_prod_prime, par.K_max- 1);
        Plotter.draw_W(W, par.K_max);

    end
    
    
    
    



    fprintf('Computing FCDS quantities.\n')
    
    tic;
    fprintf('  Construction psi function and inversion.');
    [ psi, nDS ] = build_fcds_psi(par, e_inv_des);
    fprintf('  (%0.2fs).\n', toc);
    
    tic;
    fprintf('  Computing Veta helper function.');
    [ Veta, Veta_integrand, Vexp ] = build_fcds_Veta( par, FEta );
    fprintf('  (%0.2fs).\n', toc);
    
    tic;
    fprintf('  Computing Szeta helper function.');
    [ Szeta ] = build_fcds_Szeta( par, q_prod, W);
    fprintf('  (%0.2fs).\n', toc);
    
    
    

    if (par.plot_verbosity >= 1)
        Plotter.draw_anticipation_st_curve(psi, par.T, par.K_max, par.s_c);
    end  
    if (par.plot_verbosity >= 2)
        Plotter.draw_Veta_integrand(Veta_integrand, par.K_max);
        Plotter.draw_Veta(Veta, par.K_max);
        Plotter.draw_Vexp(Vexp, par.K_max);
        Plotter.draw_Szeta(Szeta, par.K_max);
    end
    
    


    %
    % TOTAL QUIESCENT MODE
    %

    fprintf('Computing total quiescent mode pdfs.\n')
    
    tic;
    fprintf('  Computing f_qui using hvqs integral method.\n');
    [ f_qui ] = build_quiescent_mode_hvqs( par, FEta, Veta, Vexp, q_prod_prime, Szeta );
    fprintf('  Done after (%0.2fs).\n', toc);
    
    
    
    if (par.plot_verbosity >= 1)
        Plotter.draw_hvqs_quiescent_failure(f_qui, par.T);
    end


    %
    % Extreme mode
    %

    fprintf('Computing extreme mode pdfs.\n')

    tic;
    fprintf('  Computing f_ex using Poisson method.');
    [ f_ex ] = build_extreme_mode(par, fS, M_qui);
    fprintf('  (%0.2fs).\n', toc);
    
    
    
    if (par.plot_verbosity >= 2)
        Plotter.draw_extreme_failure_pdf(f_ex, par.T);
    end
    
    

    %
    % Full pdf
    %

    fprintf('Computing full pdfs.\n')

    [ fail_pdf ] = build_combined_pdf( par, f_qui, f_ex );

    if (par.plot_verbosity >= 1)
        Plotter.draw_total_failure_pdf(fail_pdf, par.T);
    end
    
    [ tail_stats ] = calculate_tail_statistics( par, f_qui, f_ex, fail_pdf );
    
    
    
    results = struct;
    results.fail_pdf = fail_pdf;
    results.f_qui = f_qui;
    results.fna = fna;
    results.q_prod = q_prod;
    results.tail_stats = tail_stats;

end

