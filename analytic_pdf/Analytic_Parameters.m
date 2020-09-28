classdef Analytic_Parameters < handle
    %ANALYTIC_PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        s_c;
        d_c;
        d_a;

        % density of plotting points

        n_s;
        n_K;   % important for final calculation!
        n_t;   % important for final calculation!
        n_d;
        n_phi;
        
        s_over_factor;
        d_over_factor;

        kappa_eps;
        sigma_eps;
        
        T_max;
        
        envelope_ascending;
        envelope_descending;
        
        load_spectrum;
        spectrum_norm_method;
        
        ds;   
        dt;
        dk;
        d_phi;
        
        K_max;
        K_min;
        delta_K;
        
        phi_max;
        
        S;
        T;
        K;
        Phi_grid;
        
      
        f_qui_abs_tol;
        f_qui_rel_tol;
        
        kappa_eta_pdf_formula;
        eta_prime_fd_algorithm;
        q_prod_prime_algorithm;
        Veta_algorithm;
        Veta_integrand_rule;
        Veta_interp_region;
        Veta_exp_rule;
        Veta_exp_clip;
        szeta_lims;
        hvqs_lims;
        
        env_des_monotonize;
        
        env_n_recovery_samples; 
        env_recover_T_rule;
        env_phi_noise_mag_rel;
        env_phi_noise_mag_abs;
        env_psi_noise_magnitude;
        env_exp_T_noise_mag_abs;
        env_exp_S_noise_mag_abs;
        
        env_des_fit_algorithm;
        env_fcds_fit_order;
        
        env_nu_samples;
        env_nu_noise_mag_abs;
        env_nu_noise_mag_rel;
        env_nu_S_noise_mag_abs;
        
        env_nu_smoothing;
        nu_envelope_class;
        nu_envelope_monotonize;
        
        nu_recovered_lin_slope;
        nu_recovered_peak_slope;
        
        verbosity;
        plot_verbosity;
    end
    
    methods
        function [ par ] = Analytic_Parameters()
            par.s_c = 1;
            par.d_c = 1;
            par.d_a = 300;

            % density of plotting points

            par.n_s = 256;
            par.n_K = 256;
            par.n_t = 128;   % important for final calculation!
            par.n_d = 128;
            par.n_phi = 128;
            
            par.s_over_factor = 1.5;
            par.d_over_factor = 100;

            par.kappa_eps = 1e-3;
            par.sigma_eps = 1e-6;
            
            %par.T_max = 150;
            par.T_max = 72000;
            
            par.K_max = 2;
            par.K_min = 1;
            
            %par.delta_K = 0.1;
            par.delta_K = (9.7e-3/250)*(par.s_c/par.d_c);
            
            par.phi_max = par.T_max;
            
            par.envelope_ascending = 'uber-ascending';
            par.envelope_descending = 'uber-descending';
           
            par.load_spectrum = 'independent-spike-model';
            
            par.spectrum_norm_method = 'integral';


            
            % for numerical integrations
            par.f_qui_abs_tol = 1e-8;
            par.f_qui_rel_tol = 1e-4;
            
            %par.eta_prime_fd_algorithm = 'method-1';
            par.eta_prime_fd_algorithm = 'method-2';
            
            par.kappa_eta_pdf_formula = 'prime-inv-combined';
            %par.kappa_eta_pdf_formula = 'prime-inv-separate';
            
            %par.q_prod_prime_algorithm = 'finite-difference';
            par.q_prod_prime_algorithm = 'analytic';
            
            %par.Veta_algorithm = 'exact';
            %par.Veta_algorithm = 'interp';
            par.Veta_algorithm = 'sum-of-antiderivatives';
            
            %par.Veta_integrand_rule = 'exact';
            par.Veta_integrand_rule = 'interp';
            
            %par.Veta_interp_region = 'zero-k-max';
            %par.Veta_interp_region = 'min-k-max';
            %par.Veta_interp_region = 'phi-k-max';
            par.Veta_interp_region = 'fat-phi-k-max';
            
            par.Veta_exp_rule = 'exponential';
            %par.Veta_exp_rule = 'taylor-expansion';
            
            par.Veta_exp_clip = true;
            
            par.szeta_lims = 'zero-k-max';
            %par.szeta_lims = 'infinity';
            
            %par.hvqs_lims = 'zero-k-max';
            %par.hvqs_lims = 'infinity';
            par.hvqs_lims = 'zero-phi';
            
            
            
            par.env_des_monotonize = true;
            
            par.env_fcds_fit_order = 3;
            
            %par.env_recover_T_rule = 'full-interval';
            par.env_recover_T_rule = 'initial-region';
            
            par.env_n_recovery_samples = 100;
            par.env_phi_noise_mag_rel = 0.05;
            par.env_phi_noise_mag_abs = 100;
            par.env_psi_noise_magnitude = 100;
            par.env_exp_T_noise_mag_abs = par.T_max/1e3;
            par.env_exp_S_noise_mag_abs = 0.02;
            
            par.env_nu_samples = 100;
            par.env_nu_noise_mag_abs = 100;
            par.env_nu_noise_mag_rel = 0.01;
            par.env_nu_S_noise_mag_abs = 0.02;
            
            par.env_nu_smoothing = 'rlowess';
            %par.env_nu_smoothing = 'none';
            
            %par.nu_envelope_class = 'linear';
            par.nu_envelope_class = 'spline';
            
            par.nu_envelope_monotonize = true;
            
            par.nu_recovered_lin_slope = 2.8;
            par.nu_recovered_peak_slope = 0;
            
            par.verbosity = 2;
            par.plot_verbosity = 1;
            
            par.update_derived_parameters();

        end
        
        function [ outcode ] = update_derived_parameters(par)
            par.phi_max = par.T_max;
            
            par.S = linspace(0, par.s_over_factor*par.s_c, par.n_s);
            par.ds = par.S(2) - par.S(1);
            par.S = par.S + par.ds/2;       % avoid some 0/0 problems when s=0 exactly

            
            par.T = linspace(1, par.T_max, par.n_t);
            par.dt = par.T(2) - par.T(1);
            
            par.K = linspace(par.K_min, par.K_max, par.n_K);
            par.dk = par.K(2) - par.K(1);
            
            par.Phi_grid = linspace(0, par.phi_max, par.n_phi);
            par.d_phi = par.Phi_grid(2) - par.Phi_grid(1);
            
            outcode = 1;
        end
        
    end
end

