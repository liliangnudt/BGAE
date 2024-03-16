function [obj] = callobj(numview, X, Up, V, A, Z, Wp, alpha, gamma, eta, lambda)
term1 = 0;
term3 = 0;

Dcol = diag(sum(Z,1));
term2 = eta * trace(V * V' - 2 * V * Z * A' + A * Dcol * A');

for p = 1 : numview
    temp1 = (X{p}-Up{p} * V)';
    temp2 = sqrt(sum(temp1.*temp1, 2));
    term1 = term1 + alpha(p)^2 * sum(temp2);
    term3 = term3 + gamma(p)^2 * norm(Wp{p} * Z' - X{p},'fro')^2;
end

obj = term1 + term2 + lambda * term3;
end