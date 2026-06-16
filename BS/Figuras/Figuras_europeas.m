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

%% COEFICIENTES AUXILIARES
R = @(t) integral(r,t,T);
Q = @(t) integral(q,t,T);

%% CALL EUROPEA
Phi_call = @(S) max(S-K,0);
g1_call  = @(t) 0;
g2_call  = @(t) S_max*exp(-Q(t)) - K*exp(-R(t));

[f_call,S_call,t_call] = bs_cn_t(0,S_max,T,N,M,r,q,sigma,Phi_call,g1_call,g2_call);

V_call_t0 = f_call(:,1);
payoff_call = Phi_call(S_call(:));

figure;
plot(S_call,V_call_t0,'LineWidth',1.5,'Color',[1 0 0]);
hold on;
plot(S_call,payoff_call,'--','LineWidth',1.2,'Color',[0 0 1]);
xlabel('S');
ylabel('f(S,0)');
grid on;
legend('Call europea','Payoff call','Location','best','FontSize', 18);
title('Call europea y payoff','FontSize', 18);

%% PUT EUROPEA
Phi_put = @(S) max(K-S,0);
g1_put  = @(t) K*exp(-R(t));
g2_put  = @(t) 0;

[f_put,S_put,t_put] = bs_cn_t(0,S_max,T,N,M,r,q,sigma,Phi_put,g1_put,g2_put);

V_put_t0 = f_put(:,1);
payoff_put = Phi_put(S_put(:));

figure;
plot(S_put,V_put_t0,'LineWidth',1.5,'Color',[1 0 0]);
hold on;
plot(S_put,payoff_put,'--','LineWidth',1.2,'Color',[0 0 1]);
xlabel('S');
ylabel('f(S,0)');
grid on;
legend('Put europea','Payoff put','Location','best','FontSize', 18);
title('Put europea y payoff','FontSize', 18);
