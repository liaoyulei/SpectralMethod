%calculate L_n'(x) and L_n(x) for the Legendre functions by
%(n+1)L_{n+1}(x)=(2n+1)xL_n(x)-nL_{n-1}(x)
%(2n+1)L_n(x)=L_{n+1}'(x)-L_{n-1}'(x)
function [dy, y] = lepoly(n, x)
if n == 0
    dy = zeros(size(x));
    y = ones(size(x));
elseif n == 1
    dy = ones(size(x));
    y = x;
else
    dy0 = zeros(size(x));
    y0 = ones(size(x));
    dyk = ones(size(x));
    yk = x;
    for k = 1: n-1
        dy = (2*k+1) * yk + dy0;
        y = ((2*k+1) * x .* yk - k * y0) / (k+1);
        dy0 = dyk;
        y0 = yk;
        dyk = dy;
        yk = y;
    end
end
end