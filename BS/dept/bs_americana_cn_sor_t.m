function [f,S,t,iteraciones,error_sor] = bs_americana_cn_sor_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2,omega,tol,max_iter)
% ==========================================================
% Valoracion de opciones americanas mediante Crank-Nicolson
% y metodo SOR proyectado
% ==========================================================
%
% ENTRADAS:
%   a,b       : extremos del dominio espacial
%   T         : vencimiento
%   N,M       : numero de subintervalos espacial y temporal
%   r,q       : tipo de interes y dividendo continuo, funciones de t
%   sigma     : volatilidad, funcion de t
%   Phi       : payoff, funcion de S
%   g1,g2     : condiciones de contorno, funciones de t
%   omega     : parametro de relajacion del metodo SOR
%   tol       : tolerancia del criterio de parada
%   max_iter  : numero maximo de iteraciones SOR por paso temporal
%
% SALIDAS:
%   f          : matriz de aproximaciones, f(i+1,j+1) ~ f(S_i,t_j)
%   S          : malla espacial
%   t          : malla temporal
%   iteraciones: iteraciones SOR realizadas en cada paso temporal
%   error_sor  : error final SOR en cada paso temporal
%
% NOTA:
%   Este metodo resuelve el problema de complementariedad asociado
%   a una opcion americana. En cada iteracion SOR se proyecta:
%
%       f_i^{k+1} = max( Phi(S_i), valor_continuacion )
%
% ==========================================================

    % ------------------------------------------------------
    % Valores por defecto
    % ------------------------------------------------------

    if nargin < 12 || isempty(omega)
        omega = 1.2;
    end

    if nargin < 13 || isempty(tol)
        tol = 1e-8;
    end

    if nargin < 14 || isempty(max_iter)
        max_iter = 10000;
    end

    % ------------------------------------------------------
    % Malla
    % ------------------------------------------------------

    h   = (b-a)/N;
    tau = T/M;

    S = linspace(a,b,N+1);
    t = linspace(0,T,M+1);

    % Matriz solucion
    f = zeros(N+1,M+1);

    % Vectores de control de convergencia
    iteraciones = zeros(M,1);
    error_sor   = zeros(M,1);

    % ------------------------------------------------------
    % Condiciones de contorno
    % ------------------------------------------------------

    for j = 1:M+1
        f(1,j)   = g1(t(j));
        f(N+1,j) = g2(t(j));
    end

    % ------------------------------------------------------
    % Condicion final
    % ------------------------------------------------------

    for i = 1:N+1
        f(i,M+1) = Phi(S(i));
    end

    % Payoff en nodos interiores
    phi = Phi(S(2:N));
    phi = phi(:);

    % ------------------------------------------------------
    % Iteracion hacia atras en el tiempo
    % ------------------------------------------------------

    for j = M:-1:1

        % --------------------------------------------------
        % Coeficientes de Crank-Nicolson
        % --------------------------------------------------

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

            % Coeficientes en t_j
            alpha_j(i) = (tau/4) * (sigmaj^2*Si^2/h^2 ...
                         - (rj-qj)*Si/h);

            beta_j(i)  = (tau/2) * (sigmaj^2*Si^2/h^2 ...
                         + rj);

            gamma_j(i) = (tau/4) * (sigmaj^2*Si^2/h^2 ...
                         + (rj-qj)*Si/h);

            % Coeficientes en t_{j+1}
            alpha_j1(i) = (tau/4) * (sigmaj1^2*Si^2/h^2 ...
                          - (rj1-qj1)*Si/h);

            beta_j1(i)  = (tau/2) * (sigmaj1^2*Si^2/h^2 ...
                          + rj1);

            gamma_j1(i) = (tau/4) * (sigmaj1^2*Si^2/h^2 ...
                          + (rj1-qj1)*Si/h);
        end

        % --------------------------------------------------
        % Construccion del segundo miembro:
        %
        %   h^{j+1} = C^{j+1} f^{j+1} + d^{j+1}
        %
        % donde el sistema de continuacion es:
        %
        %   A^j f^j = h^{j+1}
        %
        % --------------------------------------------------

        h_j1 = (1 - beta_j1) .* f(2:N,j+1);

        if N > 2
            h_j1(2:end)   = h_j1(2:end)   + alpha_j1(2:end)   .* f(2:N-1,j+1);
            h_j1(1:end-1) = h_j1(1:end-1) + gamma_j1(1:end-1) .* f(3:N,j+1);
        end

        d_j1 = zeros(N-1,1);
        d_j1(1)   = alpha_j(1)   * g1(t(j)) + alpha_j1(1)   * g1(t(j+1));
        d_j1(end) = gamma_j(end) * g2(t(j)) + gamma_j1(end) * g2(t(j+1));

        h_j1 = h_j1 + d_j1;

        % --------------------------------------------------
        % Metodo SOR proyectado
        % --------------------------------------------------
        %
        % Partimos de la solucion en la capa temporal
        % posterior f^{j+1}.
        %
        % En cada nodo:
        %
        %   \widetilde f_i^{k+1} = (1-omega) f_i^k
        %               + omega/(1+beta_i^j)
        %                 ( h_i^{j+1} + alpha_i^j f_{i-1}^{k+1}
        %                       + gamma_i^j f_{i+1}^{k} )
        %
        % y se proyecta:
        %
        %   f_i^{k+1} = max( Phi(S_i), \widetilde f_i^{k+1} )
        %
        % --------------------------------------------------

        f_k = f(2:N,j+1);
        f_k = f_k(:);

        % Proyeccion inicial sobre la restriccion americana
        f_k = max(f_k,phi);

        for k = 1:max_iter

            f_k1 = f_k;

            for i = 1:N-1

                % Nodo izquierdo ya actualizado
                if i == 1
                    termino_izq = 0;
                else
                    termino_izq = alpha_j(i) * f_k1(i-1);
                end

                % Nodo derecho aun no actualizado
                if i == N-1
                    termino_der = 0;
                else
                    termino_der = gamma_j(i) * f_k(i+1);
                end

                % Valor de continuacion aproximado por SOR
                f_tilde = (1-omega)*f_k(i) ...
                          + omega/(1+beta_j(i)) ...
                          * (h_j1(i) + termino_izq + termino_der);

                % Proyeccion sobre la restriccion americana
                f_k1(i) = max(phi(i),f_tilde);
            end

            err = norm(f_k1 - f_k,inf);

            if err < tol
                break;
            end

            f_k = f_k1;
        end

        % Guardamos informacion de convergencia
        iteraciones(j) = k;
        error_sor(j)   = err;

        if k == max_iter && err >= tol
            warning('SOR no ha convergido en t = %.6f. Error = %.3e', t(j), err);
        end

        % Solucion en la capa temporal j
        f(2:N,j) = f_k1;
    end
end
