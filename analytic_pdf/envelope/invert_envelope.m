function [ e_inv_asc, e_inv_des ] = invert_envelope( ana_par, e_forward )
%INVERT_ENVELOPE Summary of this function goes here
%   Detailed explanation goes here
    e_eps = 0.001;
    
    DD1 = linspace(0, ana_par.d_c, 100*ana_par.n_d);
    EE1 = e_forward(DD1);
    
    
    DD1 = [0, 0, DD1, ana_par.d_c, ana_par.d_c];
    EE1 = [-ana_par.s_c, -e_eps, EE1, ana_par.s_c + e_eps, 2*ana_par.s_c];
    
    e_inv_asc = @(s) interp1(EE1,DD1,s, 'linear', 'extrap');
    
    a = 4;  % switchover point
    DD2p1 = linspace(ana_par.d_c, a, 5*ana_par.n_d);
    DD2p1 = DD2p1(2:end);
    DD2p2 = logspace(log(a)/log(10), 2.5, 20*ana_par.n_d);
    DD2p2 = DD2p2(2:end);
    DD2 = [DD2p1, DD2p2];
    
    EE2 = e_forward(DD2);
    
    DD2_aug = fliplr([ana_par.d_c, ana_par.d_c, DD2]);
    EE2_aug = fliplr([2*ana_par.s_c, ana_par.s_c+e_eps, EE2]);
    
    e_inv_des = @(s) interp1(EE2_aug,DD2_aug,s, 'linear', 'extrap');
end

