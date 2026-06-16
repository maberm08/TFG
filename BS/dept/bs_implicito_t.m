function [f,S,t] = bs_implicito_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2)
% Esquema implícito con coeficientes dependientes de t
% ENTRADAS:
%   a,b    : extremos del dominio espacial
%   T      : vencimiento
%   N,M    : número de subintervalos espacial y temporal
%   r,q    : tipo de interés y dividendo continuo, funciones
%   sigma  : volatilidad, función
%   Phi    : payoff final, función
%   g1,g2  : condiciones de contorno, función
%
% SALIDAS:
%   f      : matriz de aproximaciones, f(i+1,j+1) ~ f(S_i,t_j)
%   S      : malla espacial
%   t      : malla temporal

    h   = (b-a)/N;
    tau = T/M;

    S = linspace(a,b,N+1);
    t = linspace(0,T,M+1);

    f = zeros(N+1,M+1);

    % Condiciones de contorno
    for j = 1:M+1
        f(1,j)   = g1(t(j));
        f(N+1,j) = g2(t(j));
    end

    % Condición final
    for i = 1:N+1
        f(i,M+1) = Phi(S(i));
    end

    % Iteración hacia atrás en el tiempo
    for j = M:-1:1
        alpha_j = zeros(N-1,1);
        beta_j  = zeros(N-1,1);
        gamma_j = zeros(N-1,1);

        rj     = r(t(j));
        qj     = q(t(j));
        sigmaj = sigma(t(j));

        for i = 1:N-1
            Si = S(i+1);
            alpha_j(i) = (tau/2) * (sigmaj^2*Si^2/h^2 - (rj-qj)*Si/h);
            beta_j(i)  = 1 + tau*(sigmaj^2*Si^2/h^2 + rj);
            gamma_j(i) = (tau/2) * (sigmaj^2*Si^2/h^2 + (rj-qj)*Si/h);
        end

        % Diagonales de la matriz B^j
        subdiag_B   = -alpha_j(2:end);
        diag_B      = beta_j;
        superdiag_B = -gamma_j(1:end-1);

        % Factorización LU de B^j
        [subdiag_L, diag_U, superdiag_U] = lu_tridiag(subdiag_B, diag_B, superdiag_B);

        % Término d^j de contorno
        d_j = zeros(N-1,1);
        d_j(1)   = alpha_j(1)   * g1(t(j));
        d_j(end) = gamma_j(end) * g2(t(j));

        rhs_j = f(2:N,j+1) + d_j;

        % Sustitución hacia delante: L y = rhs_j
        y = zeros(N-1,1);
        y(1) = rhs_j(1);
        for k = 2:N-1
            y(k) = rhs_j(k) - subdiag_L(k-1)*y(k-1);
        end

        % Sustitución hacia atrás: U f^j = y
        f_j = zeros(N-1,1);
        f_j(end) = y(end) / diag_U(end);
        for k = N-2:-1:1
            f_j(k) = (y(k) - superdiag_U(k)*f_j(k+1)) / diag_U(k);
        end

        f(2:N,j) = f_j;
    end
end
