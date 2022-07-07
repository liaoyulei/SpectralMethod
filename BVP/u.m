%k > 0, u = \sin(k\pi x)
%k = 0, u = \begin{cases}
%cosh(x+1) - \dfrac{x^2}2 - x, -1\le x\le 0,
%cosh(x+1) - cosh(x) - x + 1, 0\le x\le 1
%\end{cases}
function y = u(x, k)
if k > 0
    y = sin(k*pi*x);
else
    y = cosh(x+1) - x - (x < 0) .* x.^2/2 - (x > 0) .* (cosh(x) - 1);
end