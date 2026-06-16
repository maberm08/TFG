%% PUT AMERICANA Y PAYOFF
% Valora una put americana mediante Crank-Nicolson + SOR proyectado
% y representa su precio en t=0 junto con el payoff.

clear; clc; close all;

%% PARAMETROS DEL PROBLEMA

K     = 100;
T     = 1;
S_max = 300;

N = 4000;
M = 4000;

omega    = 1.2;
tol      = 1e-8;
max_iter = 10000;

%% COEFICIENTES DEPENDIENTES DEL TIEMPO

r     = @(t) 0.04 + 0.01*t;
q     = @(t) 0.01 + 0.02*t;
sigma = @(t) 0.18 + 0.04*t;

%% PAYOFF Y CONDICIONES DE CONTORNO

Phi = @(S) max(K-S,0);

% Put americana:
% En S=0 se ejerce inmediatamente y se recibe K.
% En S=S_max la put queda prácticamente sin valor.
g1_am = @(t) K;
g2_am = @(t) 0;

%% VALORACION NUMERICA

[f_am,S_am,t_am] = ...
    bs_americana_cn_sor_t(0,S_max,T,N,M,r,q,sigma,Phi,g1_am,g2_am,omega,tol,max_iter);

%% GRAFICA

figure;
plot(S_am,f_am(:,1),'LineWidth',1.5,'Color',[1 0 0]);
hold on;
plot(S_am,Phi(S_am(:)),'-.','LineWidth',1.2,'Color',[0 0 1]);

xlabel('S');
ylabel('f(S,0)');
grid on;

legend('Put americana','Payoff','Location','best');
title('Put americana con dividentos no constantes');

xlim([0 S_max]);