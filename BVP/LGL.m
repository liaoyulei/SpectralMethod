%return \{x_k,w_k\}_{k=0}^n for Legendre-Gauss-Lobatto
function [x, w] = LGL(n)
theta = (4 * (n: -1: 1)' - 1) / (4*n + 2) * pi;
sigma = (1 - (n - 1)/8/n^3 - 1/384/n^4 * (39 - 28 ./ sin(theta).^2)) .* cos(theta);
x = (sigma(1: n-1) + sigma(2: n)) ./ 2;
err = 1;
while max(abs(err)) > 1e-8
    [dy, y] = lepoly(n, x);
    err = (1 - x.^2) .* dy ./ (2*x .* dy - n*(n+1) * y);
    x = x - err;
end
x = [-1; x; 1];
[~, y] = lepoly(n, x);
w = 2/n/(n+1) ./ y.^2;
end