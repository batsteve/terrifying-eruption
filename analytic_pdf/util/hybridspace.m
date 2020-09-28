function [XX] = hybridspace(x_swap, x_max, x_n)
%HYBRIDSPACE Summary of this function goes here
%   Detailed explanation goes here

    X1 = logspace(0, log(x_swap)/log(10), ceil(x_n/2));
    X2 = linspace(x_swap, x_max, ceil(x_n/2) - 1);
    XX = [0, X1, X2(2:end)];
end

