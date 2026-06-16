clear; clc; close all;

% ==========================================================
% Valoracion de opciones americanas mediante CN + SOR proyectado
% ==========================================================
%
% Este script usa:
%   - bs_americana_cn_sor_t.m
%   - bs_cn_t.m, solo para comparar con la europea correspondiente
%
% Permite elegir:
%   - call americana
%   - put americana
%
% ==========================================================

% ----------------------------------------------------------
% 1. Seleccion del producto
% ----------------------------------------------------------

tipo_opcion = 'call';   % 'call' o 'put'

% ----------------------------------------------------------
% 2. Parametros financieros
% ----------------------------------------------------------

K  = 100;
T  = 1;
S0 = 100;

S_max = 300;

% Coeficientes dependientes del tiempo
r     = @(t) 0.05 + 0*t;
q     = @(t) 0.05 + 0*t;
sigma = @(t) 0.20 + 0*t;

% ----------------------------------------------------------
% 3. Parametros numericos
% ----------------------------------------------------------

a = 0;
b = S_max;

N = 300;
M = 400;

omega    = 1.2;
tol      = 1e-8;
max_iter = 10000;

% ----------------------------------------------------------
% 4. Payoff y condiciones de contorno americanas
% ----------------------------------------------------------

R = @(t) integral(r,t,T);
Q = @(t) integral(q,t,T);

switch tipo_opcion

    case 'call'

        Phi = @(S) max(S-K,0);

        % En S = 0, la call no tiene valor.
        g1_am = @(t) 0;

        % En S = S_max:
        g2_am = @(t) max(S_max*exp(-Q(t)) - K*exp(-R(t)), S_max - K);

    case 'put'

        Phi = @(S) max(K-S,0);

        % En S = 0, la put americana puede ejercerse inmediatamente.
        g1_am = @(t) K;

        % En S = S_max, la put queda muy out-of-the-money.
        g2_am = @(t) 0;

    otherwise

        error('tipo_opcion debe ser call o put.');
end

% ----------------------------------------------------------
% 5. Valoracion americana mediante CN + SOR proyectado
% ----------------------------------------------------------

[F_am,S_am,t_am,iteraciones,error_sor] = ...
    bs_americana_cn_sor_t(a,b,T,N,M,r,q,sigma,Phi,g1_am,g2_am,...
                          omega,tol,max_iter);

S_am_col = S_am(:);
F_am_t0  = F_am(:,1);
F_am_t0  = F_am_t0(:);

precio_am = interp1(S_am_col,F_am_t0,S0,'linear');

% ----------------------------------------------------------
% 6. Europea correspondiente para comparar
% ----------------------------------------------------------
%
% La europea NO se usa para valorar la americana.
% Solo sirve como referencia:
%   - put americana >= put europea
%   - call americana = call europea si q(t) = 0
%

switch tipo_opcion

    case 'call'

        g1_eu = @(t) 0;

        g2_eu = @(t) S_max*exp(-Q(t)) - K*exp(-R(t));

    case 'put'

        g1_eu = @(t) K*exp(-R(t));

        g2_eu = @(t) 0;
end

[F_eu,S_eu,t_eu] = bs_cn_t(a,b,T,N,M,r,q,sigma,Phi,g1_eu,g2_eu);

S_eu_col = S_eu(:);
F_eu_t0  = F_eu(:,1);
F_eu_t0  = F_eu_t0(:);

precio_eu = interp1(S_eu_col,F_eu_t0,S0,'linear');

% ----------------------------------------------------------
% 7. Prima de ejercicio anticipado
% ----------------------------------------------------------

prima_ejercicio = precio_am - precio_eu;

% ----------------------------------------------------------
% 8. Deteccion aproximada de la region de ejercicio en t = 0
% ----------------------------------------------------------

payoff_t0 = Phi(S_am_col);
ejercicio = abs(F_am_t0 - payoff_t0(:)) < 1e-6;

% ----------------------------------------------------------
% 9. Resultados
% ----------------------------------------------------------

fprintf('====================================================\n');
fprintf('Valoracion de opcion americana\n');
fprintf('Metodo: Crank-Nicolson + SOR proyectado\n');
fprintf('----------------------------------------------------\n');
fprintf('Tipo de opcion: %s americana\n', tipo_opcion);
fprintf('S0 = %.4f\n', S0);
fprintf('K  = %.4f\n', K);
fprintf('T  = %.4f\n', T);
fprintf('S_max = %.4f\n', S_max);
fprintf('N = %d, M = %d\n', N, M);
fprintf('omega = %.4f, tol = %.1e\n', omega, tol);
fprintf('----------------------------------------------------\n');
fprintf('Precio americana: %.10f\n', precio_am);
fprintf('Precio europea:   %.10f\n', precio_eu);
fprintf('Prima ejercicio:  %.10f\n', prima_ejercicio);
fprintf('----------------------------------------------------\n');
fprintf('Iteraciones SOR medias: %.2f\n', mean(iteraciones));
fprintf('Iteraciones SOR max.:   %d\n', max(iteraciones));
fprintf('Error SOR max.:         %.3e\n', max(error_sor));
fprintf('====================================================\n');

% ----------------------------------------------------------
% 10. Grafica de precios en t = 0
% ----------------------------------------------------------

figure;

plot(S_am_col,F_am_t0,'LineWidth',1.5);
hold on;

plot(S_eu_col,F_eu_t0,'--','LineWidth',1.2);

plot(S_am_col,payoff_t0,'-.','LineWidth',1.2);

xline(S0,':','S0');

xlabel('S');
ylabel('Precio');
grid on;

legend([tipo_opcion ' americana'], ...
       [tipo_opcion ' europea'], ...
       'Payoff', ...
       'S0', ...
       'Location','Best');

title(['Opcion ' tipo_opcion ' americana mediante CN + SOR proyectado']);

% ----------------------------------------------------------
% 11. Grafica de iteraciones SOR
% ----------------------------------------------------------

figure;

plot(t_am(1:end-1),iteraciones,'LineWidth',1.2);

xlabel('t');
ylabel('Iteraciones SOR');
grid on;

title('Iteraciones SOR por paso temporal');

% ----------------------------------------------------------
% 12. Region aproximada de ejercicio en t = 0
% ----------------------------------------------------------

figure;

plot(S_am_col,F_am_t0 - payoff_t0(:),'LineWidth',1.5);
hold on;

yline(0,'--');

xlabel('S');
ylabel('f(S,0) - \Phi(S)');
grid on;

title('Diferencia entre precio americano y payoff en t=0');
