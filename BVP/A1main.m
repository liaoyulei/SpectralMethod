%collocation scheme on LGL point \{x_i\}_{i=0}^N
clear;
k = 0; %k > 0 for smooth solution and 0 for limit regularity solution
err = zeros(1, 5);
for n = 1: 5
    N = 2^(n+2);
    [x, ~] = LGL(N); %N+1
    A = -LGLdiff(N)^2 + eye(N+1);
    c = [u(-1, k); zeros(N-1, 1); u(1, k)];
    b = f(x, k) - A * c;
    c(2: N) = A(2: N, 2: N) \ b(2: N);
    err(n) = max(abs(c - u(x, k)));
end
plot(x, u(x, k), x, c);