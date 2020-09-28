function [ so_par, ana_par, fS, sampler, e_forward ] = build_independent_spike_model()
%BUILD_INTEGRATOR_PARAMETERS_TEST_MODEL Summary of this function goes here
%   Detailed explanation goes here

    
    so_par = SOParameters();
    so_par.input_signal_length = 3e4;


    so_par.dc=1; so_par.sc=1; so_par.da=300;

    ana_par = Analytic_Parameters();
    ana_par.load_spectrum = 'comparison-independent-spike-model';
    [fS, fS_qui, ~, sampler ] = build_load_spectrum(ana_par);
    
    e_forward = @(d) exp(1)*ana_par.s_c/ana_par.d_c*d.*exp(-d/ana_par.d_c);
    
    [so_par.NNc,so_par.QQe] = so_par.run_counts(e_forward); 
    % note that param.QQe refers to amplitude - not the range as it is the case in the raincount alg
    
    so_par.extreme_height_threshold = 0.1;
    mu_s = integral(@(s) s.*fS_qui(s), 0, so_par.extreme_height_threshold);
    so_par.K_slope_mu = -mu_s/so_par.da;
    
    fprintf('Estimated Delta K = %0.6f.\n', so_par.K_slope_mu)
    
    

    ana_par.n_s = 256;
    ana_par.n_K = 256;
    ana_par.n_t = 256;   % important for final calculation!
    ana_par.n_d = 128;
    ana_par.n_phi = 128;
    
    
    
    ana_par.f_qui_abs_tol = 1e-16;
    ana_par.f_qui_rel_tol = 1e-13;
    
    ana_par.env_nu_samples = 1024;
    ana_par.env_nu_noise_mag_abs = 0;
    ana_par.env_nu_noise_mag_rel = 0;
    ana_par.env_nu_S_noise_mag_abs = 0;
    
    ana_par.verbosity = 2;
    ana_par.plot_verbosity = 0;

end

