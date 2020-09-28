function [f_interp] = build_interp2(f, X, Y, p)
%BUILD_INTERP2 Summary of this function goes here
%   Detailed explanation goes here

    if (nargin < 4)
        p = 0;
    end

    [XX, YY] = meshgrid(X, Y);
    ff = zeros(size(XX));
    for k = 1:size(XX, 1)
        for j = 1:size(XX, 2)
            ff(k, j) = f(XX(k, j), YY(k, j));
        end
    end
    
    f_interp = @(x, y) interp2(XX,YY,ff,x,y, 'linear', p);
end

