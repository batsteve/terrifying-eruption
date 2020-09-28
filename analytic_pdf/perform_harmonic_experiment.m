function [nu] = perform_harmonic_experiment(par, e_inv_asc, e_inv_des)
%PERFORM_HARMONIC_EXPERIMENT Summary of this function goes here
%   Detailed explanation goes here
    
    n_samples = par.env_nu_samples;
    err_std = par.env_nu_noise_mag_abs;
    err_mul = par.env_nu_noise_mag_rel;
    err_S = par.env_nu_S_noise_mag_abs;


    %QQ=[0.25:0.25:par.sc];
    SS = linspace(0, par.s_c, n_samples + 1);
    SS = SS(2:end);
    NNc = zeros(size(SS));
    
    SS_twiddle = randn(n_samples, 1)*err_S;
    
    NN_twidle =  randn(n_samples, 1).*err_std;
    NN_ftwiddle = (1 + randn(n_samples, 1)*err_mul);
    
    for i=1:length(SS)

        ddm = e_inv_asc(SS(i) + SS_twiddle(i));
        ddp = e_inv_des(SS(i) + SS_twiddle(i));
        
        %err = randn().*err_std;
        
        %NNc(i) = par.d_a*(1/ddm - 1/ddp); 
        NNc(i)= max(ceil(par.d_a*(1/ddm - 1/ddp)*NN_ftwiddle(i) + NN_twidle(i)), 1); 
    end
    
    
    switch par.env_nu_smoothing
        case 'rlowess'
            NNcs = smooth(NNc, 'rlowess');
            %NNcs = exp(smooth(log(NNc)));
            
        case 'none'
            NNcs = NNc;
    end
    
    nu = @(s) interp1(SS, NNcs, s, 'makima', 'extrap');


end

