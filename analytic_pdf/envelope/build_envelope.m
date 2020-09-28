function [e_inv_asc,e_inv_des] = build_envelope(par)
%
% Define the coherent envelope function, in terms of its two inverses
%

    d_c = par.d_c;
    s_c = par.s_c;
    switch par.envelope_ascending

        case 'uber-ascending'
            e_forward = @(d) exp(1)*par.s_c/par.d_c*d.*exp(-d/par.d_c);
            %e_raw = @(s) fminbnd(@(d) (e_forward(d) - s).^2, 0, d_c);
            e_raw = @(s) fzero(@(d) e_forward(d) - s, d_c/2);
            
            Se = linspace(0, par.s_c, par.n_s^2);
            Se = Se(2:(end-1));
            ee = zeros(1, length(Se));
            for k = 1:length(Se)
                ee(k) = e_raw(Se(k));
            end
            
            Se = [Se, par.s_c];
            ee = [ee, par.d_c];
            
            % use makima interpolation because we will eventually need
            % this derivative
            e_inv_asc = @(s) interp1(Se,ee,s, 'makima', 'extrap');

        otherwise
            warning('%s not recognized!\n');
    end

    switch par.envelope_descending            
        case 'uber-descending'
            e_forward = @(d) exp(1)*par.s_c/par.d_c*d.*exp(-d/par.d_c);
            e_raw = @(s) fminbnd(@(d) (e_forward(d) - s).^2, ...
                d_c, par.d_over_factor*d_c);
            
            Se = par.S;
            ee = zeros(length(Se), 1);
            for k = 1:length(Se)
                ee(k) = e_raw(Se(k));
            end
            e_inv_des = @(s) interp1(Se,ee,s, 'makima', 'extrap');

        otherwise
            warning('%s not recognized!\n');
    end
end

