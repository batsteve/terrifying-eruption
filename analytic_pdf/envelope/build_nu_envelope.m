function [e_inv_asc,e_inv_des] = build_nu_envelope(par, nu)
%BUILD_NU_ENVELOPE Summary of this function goes here
%   Detailed explanation goes here
    
    %
    % ascending
    %
    
    K0 = par.nu_recovered_lin_slope;
    K1 = par.nu_recovered_peak_slope;
    
    
    switch par.nu_envelope_class
        case 'linear'
    
            e_inv_asc = @(s) s/K0 + 0.001*s.^2;
            
        case 'spline'
            
            h00 = @(t) 2*t.^3 - 3*t.^2 + 1;
            h10 = @(t) t.^3 - 2*t.^2 + t;
            h01 = @(t) -2*t.^3 + 3*t.^2;
            h11 = @(t) t.^3 - t.^2;
            
            e_asc = @(d) K0*h10(d) + h01(d) + K1*h11(d);
            
            DD = linspace(0, 1, 257);
            SS = e_asc(DD);
            
            e_inv_asc = @(s) interp1(SS,DD,s, 'makima', 'extrap');
            
        otherwise
    end
    
% %     e_forward = @(d) exp(1)*par.s_c/par.d_c*d.*exp(-d/par.d_c);
% %     e_raw = @(s) fzero(@(d) e_forward(d) - s, par.d_c/2);
% %     
% %     Se = linspace(0, par.s_c, par.n_s^2);
% %     Se = Se(2:(end-1));
% %     ee = zeros(1, length(Se));
% %     for k = 1:length(Se)
% %         ee(k) = e_raw(Se(k));
% %     end
% %     
% %     Se = [Se, par.s_c];
% %     ee = [ee, par.d_c];
% %     
% %     % use makima interpolation because we will eventually need
% %     % this derivative
% %     e_inv_asc = @(s) interp1(Se,ee,s, 'makima', 'extrap');
    
    
    %
    % descending
    %
    
    Se = linspace(0, par.s_c, 2048);
    ee = zeros(length(Se), 1);
    for k = 1:length(Se)
        %ee(k) = e_raw(Se(k));
        ee(k) = 1/(e_inv_asc(Se(k))^(-1) - nu(Se(k))/(par.d_a));
    end
    
    switch par.nu_envelope_monotonize
        case true
            %ee = monotonize(ee);
            ee_m = zeros(size(ee));
            for k = 1:length(Se)
                ee_m(k) = max(ee(k:end));
            end
    
        case false
            ee_m = ee;
    end
    
    e_inv_des = @(s) interp1(Se,ee_m,s, 'makima', 'extrap');
     
end

