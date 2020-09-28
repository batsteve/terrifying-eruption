function [P_qui] = integrate_spectrum(par, fS, s_max)
%GET_P_QUIESCENT Summary of this function goes here
%   Detailed explanation goes here

    
    switch par.spectrum_norm_method
        case 'sum'
            S = par.S;
            FS = fS(S);
            P_qui = sum(FS.*par.ds);

        case 'integral'
            % bad when fS is discontinuous and spiky?
            % could be trouble at s=0 limit?  If there ever was, I don't
            % see it now.
            P_qui = integral(fS, 0, s_max);  
            
        otherwise
    end
    
end

