%% COMPARACION DE ESQUEMAS DE DIFERENCIAS FINITAS PARA BLACK-SCHOLES
% Este script compara los esquemas:
%   - explícito
%   - implícito
%   - Crank-Nicolson
% con la fórmula cerrada de Black-Scholes (blsprice)
% en el caso de una put europea vanilla.

clear; clc; close all;

%% PARAMETROS DEL PROBLEMA
K     = 100;      % strike
r     = 0.05;     % tipo de interés
q     = 0.00;     % dividendo continuo
sigma = 0.20;     % volatilidad
T     = 1.0;      % vencimiento

S0    = 100;      % precio spot donde queremos comparar
Smax  = 300;      % truncación del dominio espacial

N = 250;          % número de subintervalos espaciales 
M = 2000;         % número de subintervalos temporales

%% DEFINICION DEL PROBLEMA
% Dominio espacial [a,b]
a = 0;
b = Smax;

% Payoff terminal
Phi = @(S) max(K-S,0);

% Condiciones de contorno para put europea
g1 = @(t) K*exp(-r*(T-t));
g2 = @(t) 0;

%% PRECIO EXACTO CON BLSPRICE
[~, put_exacta] = blsprice(S0,K,r,T,sigma,q);

%% ESQUEMA EXPLICITO
[F_exp,S_exp,t_exp] = bs_explicito(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

% Tomamos el nodo más cercano a S0
[~,idx_exp] = min(abs(S_exp - S0));
put_exp = F_exp(idx_exp,1);

%% ESQUEMA IMPLICITO
[F_imp,S_imp,t_imp] = bs_implicito(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

[~,idx_imp] = min(abs(S_imp - S0));
put_imp = F_imp(idx_imp,1);

%% ESQUEMA CRANK-NICOLSON
[F_cn,S_cn,t_cn] = bs_cn(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

[~,idx_cn] = min(abs(S_cn - S0));
put_cn = F_cn(idx_cn,1);

%% ERRORES
err_exp = abs(put_exp - put_exacta);
err_imp = abs(put_imp - put_exacta);
err_cn  = abs(put_cn  - put_exacta);

%% VISUALIZACIÓN DE RESULTADOS
fprintf('---------------------------------------------\n');
fprintf('Comparacion con blsprice (put europea)\n');
fprintf('---------------------------------------------\n');
fprintf('Precio exacto Black-Scholes : %12.6f\n', put_exacta);
fprintf('\n');
fprintf('Esquema explicito          : %12.6f   Error = %.6e\n', put_exp, err_exp);
fprintf('Esquema implicito          : %12.6f   Error = %.6e\n', put_imp, err_imp);
fprintf('Esquema Crank-Nicolson     : %12.6f   Error = %.6e\n', put_cn,  err_cn);
fprintf('---------------------------------------------\n');

%% GRAFICA DE LA SOLUCION EN t=0
% No es especialmente útil pero me ha servido para visualizar mejor algunos
% resultados.
% figure;
% plot(S_exp, F_exp(:,1), 'LineWidth', 1.2); hold on;
% plot(S_imp, F_imp(:,1), 'LineWidth', 1.2);
% plot(S_cn , F_cn(:,1) , 'LineWidth', 1.2);
% 
% % Fórmula exacta de Black-Scholes sobre la malla
% put_bs_curve = zeros(size(S_cn));
% for k = 1:length(S_cn)
%     [~, put_bs_curve(k)] = blsprice(S_cn(k),K,r,T,sigma,q);
% end
% plot(S_cn, put_bs_curve, '--', 'LineWidth', 1.2);
% 
% xline(S0, ':');
% 
% legend('Explícito','Implícito','Crank-Nicolson','Black-Scholes exacta','Location','best');
% xlabel('S');
% ylabel('Precio de la put en t=0');
% title('Comparación de esquemas numéricos para Black-Scholes');
% grid on;
