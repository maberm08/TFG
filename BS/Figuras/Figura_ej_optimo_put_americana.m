clear; clc; close all;

%% PARAMETROS DEL PROBLEMA
K     = 100;
T     = 1;
S_max = 300;

r     = @(t) 0.04 + 0.01*t;
q     = @(t) 0.01 + 0.02*t;
sigma = @(t) 0.18 + 0.04*t;

N = 4000;
M = 4000;

omega    = 1.2;
tol      = 1e-8;
max_iter = 10000;
eps_ej   = 1e-8;

%% CONDICION FINAL Y CONTORNO
Phi   = @(S) max(K-S,0);
g1_am = @(t) K;
g2_am = @(t) 0;

[f_am,S_am,t_am,iteraciones,error_sor] = ...
    bs_americana_cn_sor_t(0,S_max,T,N,M,r,q,sigma,Phi,g1_am,g2_am,omega,tol,max_iter);

payoff = Phi(S_am(:));
S_ej_opt = nan(size(t_am));

for j = 1:length(t_am)
    idx = find(f_am(:,j) > payoff + eps_ej,1,'first');
    if ~isempty(idx)
        S_ej_opt(j) = S_am(idx);
    end
end

figure;
plot(t_am,S_ej_opt,'LineWidth',1.5,'Color',[1 0 0]);
xlabel('t');
ylabel('S');
grid on;
title('Curva de ejercicio optimo de la put americana');
