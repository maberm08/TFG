function [F,S,t] = bs_cn(a,b,T,N,M,r,q,sigma,Phi,g1,g2)
% ENTRADAS:
%   a,b    : extremos del dominio espacial
%   T      : vencimiento
%   N,M    : número de subintervalos espaciales y temporales
%   r,q    : tipo de interés y dividendo continuo (constantes)
%   sigma  : volatilidad
%   Phi    : payoff terminal, función handle
%   g1,g2  : condiciones de contorno, funciones handle
%
% SALIDAS:
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

    % Diagonales de la matriz I+B
    subdiag_IB   = -alpha(2:end);
    diag_IB      = 1 + beta;
    superdiag_IB = -gamma(1:end-1);

    % Factorización LU de I+B
    [subdiag_L, diag_U, superdiag_U] = lu_tridiag(subdiag_IB, diag_IB, superdiag_IB);

    % Iteración hacia atrás en el tiempo
    for j = M:-1:1
        % Término 2F^{j+1}
        segundo_miembro = 2 * F(2:N,j+1);

        % Término d^{j,j+1} de contorno
        d_jj1 = zeros(N-1,1);
        d_jj1(1)   = alpha(1)   * (g1(t(j)) + g1(t(j+1)));
        d_jj1(end) = gamma(end) * (g2(t(j)) + g2(t(j+1)));

        segundo_miembro = segundo_miembro + d_jj1;

        % Sustitución hacia delante: L y = segundo_miembro
        y = zeros(N-1,1);
        y(1) = segundo_miembro(1);
        for k = 2:N-1
            y(k) = segundo_miembro(k) - subdiag_L(k-1)*y(k-1);
        end

        % Sustitución hacia atrás: U V = y
        V = zeros(N-1,1);
        V(end) = y(end) / diag_U(end);
        for k = N-2:-1:1
            V(k) = (y(k) - superdiag_U(k)*V(k+1)) / diag_U(k);
        end

        % Recuperar F^j a partir de V = F^j + F^{j+1}
        F(2:N,j) = V - F(2:N,j+1);
    end
end
