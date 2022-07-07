%spectral-Galerkin scheme scheme on LGL point \{x_i\}_{i=0}^N
clear;
k = 0; %k > 0 for smooth solution and 0 for limit regularity solution
err = zeros(1, 5);
for n = 1: 5
    N = 2^(n+2);
    [x, w] = LGL(N); %N+1
    gamma = [2./(1: 2: 2*N-1)'; 2/N]; %N+1
    A = diag(1 + 2./(1: 2: 2*N-3)./(5: 2: 2*N+1)) - diag(1./sqrt((3: 2: 2*N-5).*(7: 2: 2*N-1))./(5: 2: 2*N-3), -2) - diag(1./sqrt((3: 2: 2*N-5).*(7: 2: 2*N-1))./(5: 2: 2*N-3), 2);
    Q = zeros(N+1); %Q_{ij}=L_j(x_i)
    Q(:, [1, 2]) = [ones(N+1, 1), x];
    for j = 1: N-1
        Q(:, j+2) = ((2*j+1) * x .* Q(:, j+1) - j * Q(:, j)) / (j+1);
    end
    b = (Q'./gamma) * (w.*f(x, k)); %\hat{f}=diag(1/gamma)Q'diag(w)f
    b = b .* gamma;
    b = 1./sqrt(6: 4: 4*N-2)' .* (b(3: N+1) - b(1: N-1));
    alpha = (u(1, k) + u(-1, k)) / 2; %boundary
    beta = (u(1, k) - u(-1, k)) / 2;
    b(1) = b(1) + 2/sqrt(6) * alpha;
    b(2) = b(2) + 2/3/sqrt(10) * beta;
    c = A \ b; 
    c = [alpha-c(1)/sqrt(6); beta-c(2)/sqrt(10); c(1: N-3)./sqrt(6: 4: 4*N-10)' - c(3: N-1)./sqrt(14: 4: 4*N-2)'; c(N-2)/sqrt(4*N-6); c(N-1)/sqrt(4*N-2)];
    c = Q * c; %u=Q\hat{u}
    err(n) = max(abs(c - u(x, k)));
end
plot(x, u(x, k), x, c);