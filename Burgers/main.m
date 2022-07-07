%spectral-Galerkin scheme scheme on LGL point \{x_i\}_{i=0}^N
%solve \partial_tu+u\partial_xu=\epsilon\partial_x^2u  x\in[-1,1]
%u(\pm 1,t)=0 t\ge 0               u(x,0)=-\sin(\pi x) x\in[-1,1]
clear;
epsilon = 0.02;
T = 0.995;
t = 1e-3;
N = 128;
[x, w, Q] = LGL(N); %N+1
gamma = [2./(1: 2: 2*N-1)'; 2/N]; %N+1
A = diag(2./(1: 2: 2*N-3)./(5: 2: 2*N+1)) - diag(1./sqrt((3: 2: 2*N-5).*(7: 2: 2*N-1))./(5: 2: 2*N-3), -2) - diag(1./sqrt((3: 2: 2*N-5).*(7: 2: 2*N-1))./(5: 2: 2*N-3), 2);
A = epsilon * t * eye(N - 1) + A;
let = @(f) (Q'./gamma) * (w.*f); %\hat{f}=diag(1/gamma)Q'diag(w)f
u = zeros(N+1, T/t+1);
u(:, 1) = -sin(pi*x);
c = let(u(:, 1));
for j = 1: T/t
    b = c .* gamma;
    c = let(u(:, j).^2);
    b = 1./sqrt(6: 4: 4*N-2)' .* (b(3: N+1) - b(1: N-1)) + t * c(2: N) ./ sqrt(6: 4: 4*N-2)';
    c = A \ b;
    c = [-c(1)/sqrt(6); -c(2)/sqrt(10); c(1: N-3)./sqrt(6: 4: 4*N-10)' - c(3: N-1)./sqrt(14: 4: 4*N-2)'; c(N-2)/sqrt(4*N-6); c(N-1)/sqrt(4*N-2)];
    u(:, j+1) = Q * c;    
end
subplot(1, 2, 1);
[X, Y] = meshgrid(x, 0: t: T);
surf(X, Y, u');
subplot(1, 2, 2);
plot(x, u(:, T/t+1));