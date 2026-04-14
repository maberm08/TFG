function [F,S,t] = bs_implicito(a,b,T,N,M,r,q,sigma,Phi,g1,g2)
% ENTRADAS:
%   a,b    : extremos del dominio espacial
%   T      : vencimiento
%   N,M    : número de subintervalos espaciales y temporales
%   r,q    : tipo de interés y dividendo continuo (constantes)
%   sigma  : volatilidad
%   Phi    : payoff terminal, función
%   g1,g2  : condiciones de contorno, funciones
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

    % Coeficientes del esquema implícito
    alpha = zeros(N-1,1);
    beta  = zeros(N-1,1);
    gamma = zeros(N-1,1);

    for i = 1:N-1
        Si = S(i+1);
        alpha(i) = (tau/2) * (sigma^2*Si^2/h^2 - (r-q)*Si/h);
        beta(i)  = 1 + tau*(sigma^2*Si^2/h^2 + r);
        gamma(i) = (tau/2) * (sigma^2*Si^2/h^2 + (r-q)*Si/h);
    end

    % Diagonales de la matriz B
    subdiag_B = -alpha(2:end);
    diag_B    = beta;
    superdiag_B = -gamma(1:end-1);

    % Factorización LU de B
    [subdiag_L, diag_U, superdiag_U] = lu_tridiag(subdiag_B, diag_B, superdiag_B);

    % Iteración hacia atrás en el tiempo
    for j = M:-1:1
        segundo_miembro = F(2:N,j+1);

        % Término d^j de contorno
        d_j = zeros(N-1,1);
        d_j(1)   = alpha(1)   * g1(t(j));
        d_j(end) = gamma(end) * g2(t(j));

        segundo_miembro = segundo_miembro + d_j;

        % Sustitución hacia delante: L y = segundo_miembro
        y = zeros(N-1,1);
        y(1) = segundo_miembro(1);
        for k = 2:N-1
            y(k) = segundo_miembro(k) - subdiag_L(k-1)*y(k-1);
        end

        % Sustitución hacia atrás: U solucion_interior = y
        solucion_interior = zeros(N-1,1);
        solucion_interior(end) = y(end) / diag_U(end);
        for k = N-2:-1:1
            solucion_interior(k) = (y(k) - superdiag_U(k)*solucion_interior(k+1)) / diag_U(k);
        end

        F(2:N,j) = solucion_interior;
    end
end
