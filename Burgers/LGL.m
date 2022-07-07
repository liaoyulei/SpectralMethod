%return \{x_k,w_k\}_{k=0}^n for Legendre-Gauss-Lobatto
%Q_{ij}=L_j(x_i) for i,j=0,\cdots,n
function [x, w, Q] = LGL(n)
theta = (4 * (n: -1: 1)' - 1) / (4*n + 2) * pi;
sigma = (1 - (n - 1)/8/n^3 - 1/384/n^4 * (39 - 28 ./ sin(theta).^2)) .* cos(theta);
x = (sigma(1: n-1) + sigma(2: n)) ./ 2;
err = 1;
while max(abs(err)) > 1e-8
    [dQ, Q] = lepoly(n, x);
    err = (1 - x.^2) .* dQ(:, n+1) ./ (2*x .* dQ(:, n+1) - n*(n+1) * Q(:, n+1));
    x = x - err;
end
x = [-1; x; 1];
[~, Q] = lepoly(n, x);
w = 2/n/(n+1) ./ Q(:, n+1).^2;
end