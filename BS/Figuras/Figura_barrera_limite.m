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

num_barreras = 5;

%% COEFICIENTES AUXILIARES
R = @(t) integral(r,t,T);

%% PUT VANILLA
Phi_vanilla = @(S) max(K-S,0);
g1_vanilla  = @(t) K*exp(-R(t));
g2_vanilla  = @(t) 0;

[f_vanilla,S_vanilla,t_vanilla] = bs_cn_t(0,S_max,T,N,M,r,q,sigma,Phi_vanilla,g1_vanilla,g2_vanilla);

V_vanilla_t0 = f_vanilla(:,1);
payoff_vanilla = Phi_vanilla(S_vanilla(:));

%% DOWN-AND-OUT PUTS
S_b_vals = K ./ 2.^(1:num_barreras);
colores = lines(num_barreras);
leyendas = cell(num_barreras + 2,1);

figure;
plot(S_vanilla,V_vanilla_t0,'k','LineWidth',1.8);
hold on;
leyendas{1} = 'Put vanilla';

for n = 1:num_barreras
    S_b = S_b_vals(n);

    Phi_dop = @(S) (S > S_b) .* max(K-S,0);
    g1_dop  = @(t) 0;
    g2_dop  = @(t) 0;

    [f_dop,S_dop,t_dop] = bs_cn_t(S_b,S_max,T,N,M,r,q,sigma,Phi_dop,g1_dop,g2_dop);

    V_dop_t0 = f_dop(:,1);

    plot(S_dop,V_dop_t0,'Color',colores(n,:),'LineWidth',1.1);
    leyendas{n+1} = sprintf('DOP, S_b = %.5g', S_b);
end

plot(S_vanilla,payoff_vanilla,'k--','LineWidth',1.2);
leyendas{end} = 'Payoff put vanilla';

xlabel('S');
ylabel('f(S,0)');
grid on;
legend(leyendas,'Location','best');
title('Down-and-out puts cuando S_b tiende a 0');
