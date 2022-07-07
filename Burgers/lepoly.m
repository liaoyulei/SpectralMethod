%calculate L_j'(x_i) and L_j(x_i) j=0,\dots n for the Legendre functions by
%(n+1)L_{n+1}(x)=(2n+1)xL_n(x)-nL_{n-1}(x)
%(2n+1)L_n(x)=L_{n+1}'(x)-L_{n-1}'(x)
function [dQ, Q] = lepoly(n, x)
dQ = zeros(size(x, 1), n+1);
Q = zeros(size(x, 1), n+1);
Q(:, 1) = 1;
if n > 0
    dQ(:, 2) = 1;
    Q(:, 2) = x;
    for j = 1: n-1
        dQ(:, j+2) = (2*j+1) * Q(:, j+1) + dQ(:, j);
        Q(:, j+2) = ((2*j+1) * x .* Q(:, j+1) - j * Q(:, j)) / (j+1);
    end
end
end