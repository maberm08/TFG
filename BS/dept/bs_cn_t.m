function [f,S,t] = bs_cn_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2)
% Esquema de Crank-Nicolson con coeficientes dependientes de t
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
        alpha_j  = zeros(N-1,1);
        beta_j   = zeros(N-1,1);
        gamma_j  = zeros(N-1,1);

        alpha_j1 = zeros(N-1,1);
        beta_j1  = zeros(N-1,1);
        gamma_j1 = zeros(N-1,1);

        rj      = r(t(j));
        qj      = q(t(j));
        sigmaj  = sigma(t(j));

        rj1     = r(t(j+1));
        qj1     = q(t(j+1));
        sigmaj1 = sigma(t(j+1));

        for i = 1:N-1
            Si = S(i+1);

            % Coeficientes en t_j  (parte implícita)
            alpha_j(i) = (tau/4) * (sigmaj^2*Si^2/h^2 - (rj-qj)*Si/h);
            beta_j(i)  = (tau/2) * (sigmaj^2*Si^2/h^2 + rj);
            gamma_j(i) = (tau/4) * (sigmaj^2*Si^2/h^2 + (rj-qj)*Si/h);

            % Coeficientes en t_{j+1} (parte explícita)
            alpha_j1(i) = (tau/4) * (sigmaj1^2*Si^2/h^2 - (rj1-qj1)*Si/h);
            beta_j1(i)  = (tau/2) * (sigmaj1^2*Si^2/h^2 + rj1);
            gamma_j1(i) = (tau/4) * (sigmaj1^2*Si^2/h^2 + (rj1-qj1)*Si/h);
        end

        % Diagonales de A^j
        subdiag_A   = -alpha_j(2:end);
        diag_A      = 1 + beta_j;
        superdiag_A = -gamma_j(1:end-1);

        % Factorización LU de A^j
        [subdiag_L, diag_U, superdiag_U] = lu_tridiag(subdiag_A, diag_A, superdiag_A);

        h_j1 = (1 - beta_j1) .* f(2:N,j+1);

        if N > 2
            h_j1(2:end)   = h_j1(2:end)   + alpha_j1(2:end)   .* f(2:N-1,j+1);
            h_j1(1:end-1) = h_j1(1:end-1) + gamma_j1(1:end-1) .* f(3:N,j+1);
        end

        d_j1 = zeros(N-1,1);
        d_j1(1)   = alpha_j(1)   * g1(t(j))   + alpha_j1(1)   * g1(t(j+1));
        d_j1(end) = gamma_j(end) * g2(t(j))   + gamma_j1(end) * g2(t(j+1));

        h_j1 = h_j1 + d_j1;

        % Sustitución hacia delante: L y = h_j1
        y = zeros(N-1,1);
        y(1) = h_j1(1);
        for k = 2:N-1
            y(k) = h_j1(k) - subdiag_L(k-1)*y(k-1);
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
