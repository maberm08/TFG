function [F,S,t] = bs_implicito_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2)
% Esquema explícito con coeficientes dependientes de t
% ENTRADAS:
%   a,b    : extremos del dominio espacial
%   T      : vencimiento
%   N,M    : número de subintervalos espacial y temporal
%   r,q    : tipo de interés y dividendo continuo, funciones
%   sigma  : volatilidad, función
%   Phi    : payoff terminal, función
%   g1,g2  : condiciones de contorno, función
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

    % Iteración hacia atrás en el tiempo
    for j = M:-1:1
        alpha = zeros(N-1,1);
        beta  = zeros(N-1,1);
        gamma = zeros(N-1,1);

        rj     = r(t(j));
        qj     = q(t(j));
        sigmaj = sigma(t(j));

        for i = 1:N-1
            Si = S(i+1);
            alpha(i) = (tau/2) * (sigmaj^2*Si^2/h^2 - (rj-qj)*Si/h);
            beta(i)  = 1 + tau*(sigmaj^2*Si^2/h^2 + rj);
            gamma(i) = (tau/2) * (sigmaj^2*Si^2/h^2 + (rj-qj)*Si/h);
        end

        % Diagonales de la matriz B^j
        subdiag_B   = -alpha(2:end);
        diag_B      = beta;
        superdiag_B = -gamma(1:end-1);

        % Factorización LU de B^j
        [subdiag_L, diag_U, superdiag_U] = lu_tridiag(subdiag_B, diag_B, superdiag_B);

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