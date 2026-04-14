function [F,S,t] = bs_explicito(a,b,T,N,M,r,q,sigma,Phi,g1,g2)
% ENTRADAS:
%   a,b    : extremos del dominio espacial
%   T      : vencimiento
%   N,M    : número de subintervalos espacial y temporal
%   r,q    : tipo de interés y dividendo continuo (ctes)
%   sigma  : volatilidad
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

    % Matriz solución
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

    % Coeficientes del esquema explícito
    aa = zeros(N-1,1);
    bb = zeros(N-1,1);
    cc = zeros(N-1,1);

    for i = 1:N-1
        Si = S(i+1);  % corresponde a S_i, i=1,...,N-1
        aa(i) = (tau/2) * (sigma^2*Si^2/h^2 - (r-q)*Si/h);
        bb(i) = 1 - tau*(sigma^2*Si^2/h^2 + r);
        cc(i) = (tau/2) * (sigma^2*Si^2/h^2 + (r-q)*Si/h);
    end

    % Iteración hacia atrás en el tiempo
    for j = M:-1:1
        F(2:N,j) = aa .* F(1:N-1,j+1) + bb .* F(2:N,j+1) + cc .* F(3:N+1,j+1);
    end
end
