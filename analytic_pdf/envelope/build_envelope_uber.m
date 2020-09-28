function [e_inv_asc,e_inv_des] = build_envelope_uber( ana_par )
%BUILD_ENVELOPE_UBER Summary of this function goes here
%   Detailed explanation goes here

    d_eps = 0.001;
    e_eps = 0.001;
    
    e_forward = @(d) exp(1)*ana_par.s_c/ana_par.d_c*d.*exp(-d/ana_par.d_c);
    
    DD1 = linspace(0, ana_par.d_c, 100*ana_par.n_d);
    EE1 = e_forward(DD1);
    
    %DD1 = [-1, -d_eps, DD1, par.d_c+d_eps, 2*par.d_c];
    %EE1 = [0, 0, EE1, 1, 1];
    
    DD1 = [0, 0, DD1, ana_par.d_c, ana_par.d_c];
    EE1 = [-ana_par.s_c, -e_eps, EE1, ana_par.s_c + e_eps, 2*ana_par.s_c];
    
    e_inv_asc = @(s) interp1(EE1,DD1,s, 'linear', 'extrap');
    
    
    
    %DD2 = linspace(ana_par.d_c, 20*ana_par.d_c, 100*ana_par.n_d);
    
    
    a = 4;  % switchover point
    DD2p1 = linspace(ana_par.d_c, a, 5*ana_par.n_d);
    DD2p1 = DD2p1(2:end);
    DD2p2 = logspace(log(a)/log(10), 2.5, 20*ana_par.n_d);
    DD2p2 = DD2p2(2:end);
    DD2 = [DD2p1, DD2p2];
    
    EE2 = e_forward(DD2);
    
    DD2_aug = fliplr([1, 1, DD2]);
    EE2_aug = fliplr([2*ana_par.s_c, ana_par.s_c+e_eps, EE2]);
    
    e_inv_des = @(s) interp1(EE2_aug,DD2_aug,s, 'linear', 'extrap');
end

