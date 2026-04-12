function [F,S,t] = bs_cn(a,b,T,N,M,r,q,sigma,Phi,g1,g2)
% ENTRADAS:
%   a,b    : extremos del dominio espacial
%   T      : vencimiento
%   N,M    : número de subintervalos espacial y temporal
%   r,q    : tipo de interés y dividendo continuo (ctes)
%   sigma  : volatilidad
%   Phi    : payoff terminal, función handle
%   g1,g2  : condiciones de contorno, función handle
%
% SALIDA:
%   F      : matriz de aproximaciones, F(i+1,j+1) ~ f(S_i,t_j)
%   S      : malla espacial
%   t      : malla temporal

    h   = (b-a)/N;
    tau = T/M;

    S = linspace(a,b,N+1);
    t = linspace(0,T,M+1);

    F = zeros(N+1,M+1);

    % Condiciones de contorno
    for j = 1:M+1
        F(1,j)   = g1(t(j));
        F(N+1,j) = g2(t(j));
    end

    % Condición terminal
    for i = 1:N+1
        F(i,M+1) = Phi(S(i));
    end

    % Coeficientes del esquema de Crank-Nicolson
    alpha = zeros(N-1,1);
    beta  = zeros(N-1,1);
    gamma = zeros(N-1,1);

    for i = 1:N-1
        Si = S(i+1);
        alpha(i) = (tau/4) * (sigma^2*Si^2/h^2 - (r-q)*Si/h);
        beta(i)  = (tau/2) * (sigma^2*Si^2/h^2 + r);
        gamma(i) = (tau/4) * (sigma^2*Si^2/h^2 + (r-q)*Si/h);
    end

    % Matriz B de la notación A=I+B, C=I-B
    B = diag(beta) + diag(-gamma(1:N-2),1) + diag(-alpha(2:N-1),-1);

    I = eye(N-1);
    A = I + B;
    C = I - B;

    % Iteración hacia atrás en el tiempo
    for j = M:-1:1
        rhs = C * F(2:N,j+1);

        % Término de contorno b^{j,j+1}
        rhs(1)   = rhs(1)   + alpha(1)   * (g1(t(j)) + g1(t(j+1)));
        rhs(end) = rhs(end) + gamma(end) * (g2(t(j)) + g2(t(j+1)));

        % Resolver sistema
        F(2:N,j) = A \ rhs;
    end
end