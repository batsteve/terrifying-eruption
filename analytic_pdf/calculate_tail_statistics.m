function [ tail_stats ] = calculate_tail_statistics( ana_par, f_qui, f_ex, f_nfail )
%CALCULATE_TAIL_STATISTICS Summary of this function goes here
%   Detailed explanation goes here


% alef   -- mean t_fail
% bet    -- variance t_fail
% gimmel -- long tail plateau value
% dalet  -- crossover value
% he     -- near tail slope

    TT = linspace(0, ana_par.T_max, 1024);
    
    alef = integral(@(t) t.*f_nfail(t), 0, ana_par.T_max);
    
    bet = integral(@(t) (t - alef).^2.*f_nfail(t), 0, ana_par.T_max);
    
    gimmel = mean(f_ex(TT));
    
    dalet = fzero(@(t) f_qui(t) - gimmel, (1/3)*ana_par.T_max);
    
    XX = linspace(dalet, dalet + ana_par.T_max/3, 21);
    XX = XX(2:end);
    ZZ = f_qui(XX);
    EE = log((ZZ - gimmel)./gimmel)./(XX - dalet);
    
    he = median(EE);
    
    tail_stats = struct;
    tail_stats.alef = alef;
    tail_stats.bet = bet;
    tail_stats.gimmel = gimmel;
    tail_stats.dalet = dalet;
    tail_stats.he = he;
    
    
end

