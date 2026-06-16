function [f,S,t] = bs_explicito_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2)
% Esquema explícito con coeficientes dependientes de t
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

    % Matriz solución
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
        a_j1 = zeros(N-1,1);
        b_j1 = zeros(N-1,1);
        c_j1 = zeros(N-1,1);
    
        rj     = r(t(j+1));
        qj     = q(t(j+1));
        sigmaj = sigma(t(j+1));
    
        for i = 1:N-1
            Si = S(i+1);  % corresponde a S_i, i=1,...,N-1
            a_j1(i) = (tau/2) * (sigmaj^2*Si^2/h^2 - (rj-qj)*Si/h);
            b_j1(i) = 1 - tau*(sigmaj^2*Si^2/h^2 + rj);
            c_j1(i) = (tau/2) * (sigmaj^2*Si^2/h^2 + (rj-qj)*Si/h);
        end
    
        f(2:N,j) = a_j1 .* f(1:N-1,j+1) + b_j1 .* f(2:N,j+1) + c_j1 .* f(3:N+1,j+1);
    end
end
