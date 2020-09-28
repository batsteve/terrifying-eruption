function [N] = run_rainflow(param, DS)
%RUN_RAINFLOW_SJ Summary of this function goes here
%   Detailed explanation goes here

    ss = [DS.sspp(:)'; DS.ssnp(:)'];
    ss = ss(:);
    
    %
    % bisection search!
    %
    
    g1 = 1;
    g2 = length(DS.sspp(:));
    
    if (f_rainflow(ss(1:(2*g2))) < 1)
        N = g2;
        
        warning('Rainflow failure never occurs!')
    else
        while (g2 - g1) > 3
            g3 = ceil((g2 + g1)/2);
            C=f_rainflow(ss(1:(2*g3)));
            if (C > 1)
                g2 = g3;
            else
                g1 = g3;
            end
        end
        N = ceil((g2 + g1)/2);    
    end
    
    function [C] = f_rainflow(ss)
        [c,hist,edges,rmm,idx] = rainflow(ss); % here edge is range! not amplitude
        NNci = interp1(param.QQe,param.NNc,edges(2:end),'linear','extrap'); % we devide edges with 2 to transform to amplitude
        NNr=sum(hist,2);
        NNc = NNr.*NNci.^-1;
        C=sum(NNc);
    end
    
    
end

