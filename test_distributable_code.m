addpath('analytic_pdf');
addpath('analytic_pdf/mdas_methods');
addpath('analytic_pdf/fcds_methods');
addpath('analytic_pdf/util');
addpath('analytic_pdf/envelope');
addpath('time_series_integrator');


[ so_par, ana_par, fS, sampler, e_forward ] = build_independent_spike_model();


ana_par.delta_K = estimate_delta_k(ana_par, fS);
ana_par.T_max = (2.7/ana_par.delta_K) * 1.2;  % overestimation is okay
ana_par.update_derived_parameters();

ana_par.verbosity = 2;
ana_par.plot_verbosity = 1;

[ e_inv_asc, e_inv_des ] = invert_envelope( ana_par, e_forward );
[kappa, kappa_inv, kappa_prime, kappa_prime_inv] = build_kappas(ana_par, e_inv_asc);

[ ana_results ] = build_full_pdf(ana_par, fS, e_forward);


seed = nan;  % default to random seed
    
N = 10;
NN1 = zeros(N, 1);
NN2 = zeros(N, 1);

for k = 1:N
    [DS_run] = build_independent_spike_signal(so_par, seed, sampler);
    [n_full, eval_struct_full] = so_integrate_full(so_par, DS_run, e_forward, kappa);
    NN1(k) = n_full;
    [n_fast, eval_struct_full] = so_decomp_synth(so_par, DS_run, e_forward, kappa);
    NN2(k) = n_fast;
end

figure(101);
clf;
histogram(NN1);

figure(102);
clf;
histogram(NN2);

figure(103);
clf;
plot(ana_par.T, ana_results.fail_pdf(ana_par.T));
set(gca, 'YScale', 'log');
ylim([1e-8, 1e-2])