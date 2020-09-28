function [f_mono] = monotonize(f, range)
%MONOTONIZE Summary of this function goes here
%   Detailed explanation goes here

    X = linspace(range(1), range(2), 256);
    F = f(X);
    %D = diff(X);
    F_m = zeros(length(F), 1);
    
    for k = 1:length(F)
        F_m(k) = max(F(k:end));
    end
    
    f_mono = @(x) interp1(X, F_m, x, 'linear', 'extrap');

end

