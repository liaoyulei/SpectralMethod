%k > 0, f = (k^2\pi^2+1)\sin(k\pi x)
%k = 0, f = \begin{cases}
%1 - x - \dfrac{x^2}2, -1\le x\le 0,
%1 - x,                 0\le x\le 1
%\end{cases}
function y = f(x, k)
if k
    y = (k^2*pi^2 + 1) * sin(k*pi*x);
else
    y = 1 - x - (x < 0) .* x.^2/2;
end