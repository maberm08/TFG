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

%% CASO 1: SIN DIVIDENDOS
r     = @(t) 0.04 + 0.01*t;
q     = @(t) 0.00 + 0*t;
sigma = @(t) 0.18 + 0.04*t;

R = @(t) integral(r,t,T);
Q = @(t) integral(q,t,T);

Phi   = @(S) max(S-K,0);
g1_am = @(t) 0;
g2_am = @(t) max(S_max*exp(-Q(t)) - K*exp(-R(t)), S_max - K);

g1_eu = @(t) 0;
g2_eu = @(t) S_max*exp(-Q(t)) - K*exp(-R(t));

[f_am_0,S_am_0,t_am_0] = ...
    bs_americana_cn_sor_t(0,S_max,T,N,M,r,q,sigma,Phi,g1_am,g2_am,omega,tol,max_iter);

[f_eu_0,S_eu_0,t_eu_0] = ...
    bs_cn_t(0,S_max,T,N,M,r,q,sigma,Phi,g1_eu,g2_eu);

figure;
plot(S_am_0,f_am_0(:,1),'LineWidth',1.5,'Color',[1 0 0]);
hold on;
plot(S_eu_0,f_eu_0(:,1),'--','LineWidth',1.5,'Color',[0.25 0.25 0.25]);
plot(S_am_0,Phi(S_am_0(:)),'-.','LineWidth',1.2,'Color',[0 0 1]);
xlabel('S');
ylabel('f(S,0)');
grid on;
legend('Call americana','Call europea','Payoff','Location','best');
title('Call americana y europea sin dividendos');

%% CASO 2: CON DIVIDENDOS NO CONSTANTES
r     = @(t) 0.04 + 0.01*t;
q     = @(t) 0.01 + 0.02*t;
sigma = @(t) 0.18 + 0.04*t;

R = @(t) integral(r,t,T);
Q = @(t) integral(q,t,T);

g1_am = @(t) 0;
g2_am = @(t) max(S_max*exp(-Q(t)) - K*exp(-R(t)), S_max - K);

g1_eu = @(t) 0;
g2_eu = @(t) S_max*exp(-Q(t)) - K*exp(-R(t));

[f_am_q,S_am_q,t_am_q] = ...
    bs_americana_cn_sor_t(0,S_max,T,N,M,r,q,sigma,Phi,g1_am,g2_am,omega,tol,max_iter);

[f_eu_q,S_eu_q,t_eu_q] = ...
    bs_cn_t(0,S_max,T,N,M,r,q,sigma,Phi,g1_eu,g2_eu);

figure;
plot(S_am_q,f_am_q(:,1),'LineWidth',1.5,'Color',[1 0 0]);
hold on;
plot(S_am_q,Phi(S_am_q(:)),'-.','LineWidth',1.2,'Color',[0 0 1]);
xlabel('S');
ylabel('f(S,0)');
grid on;
legend('Call americana','Payoff','Location','best');
title('Call americana con dividendos no constantes');