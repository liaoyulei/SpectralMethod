%%D_{ij}=h_j'(x_i) for LGL point
function D = LGLdiff(n)
[x, ~] = LGL(n); %n+1
[~, y] = lepoly(n, x);
D = y ./ y' ./ (x - x' + eye(n+1)) - eye(n+1); %D_{ij}=\dfrac{y_i}{y_j(x_i-x_y)}, i\neq j, 0,i=j
D(1, 1) = -n*(n+1)/4;
D(n+1, n+1) = n*(n+1)/4;
end