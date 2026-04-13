%% COMPARACION DE ESQUEMAS DE DIFERENCIAS FINITAS PARA BLACK-SCHOLES
% Este script compara los esquemas:
%   - explícito
%   - implícito
%   - Crank-Nicolson
% con la fórmula cerrada de Black-Scholes (blsprice)
% en el caso de una call europea vanilla.

clear; clc; close all;

%% PARAMETROS DEL PROBLEMA
K     = 100;      % strike
r     = 0.05;     % tipo de interés
q     = 0.00;     % dividendo continuo
sigma = 0.20;     % volatilidad
T     = 1.0;      % vencimiento

S0    = 100;      % precio spot donde queremos comparar
Smax  = 300;      % truncación del dominio espacial

N = 300;          % número de subintervalos espaciales 
M = 4000;         % número de subintervalos temporales

%% DEFINICION DEL PROBLEMA
% Dominio espacial [a,b]
a = 0;
b = Smax;

% Payoff terminal
Phi = @(S) max(S-K,0);

% Condiciones de contorno para call europea
g1 = @(t) 0;
g2 = @(t) Smax*exp(-q*(T-t)) - K*exp(-r*(T-t));

%% PRECIO EXACTO CON BLSPRICE
[call_exacta, ~] = blsprice(S0,K,r,T,sigma,q);

%% ESQUEMA EXPLICITO
[F_exp,S_exp,t_exp] = bs_explicito(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

% Tomamos el nodo más cercano a S0
[~,idx_exp] = min(abs(S_exp - S0));
call_exp = F_exp(idx_exp,1);

%% ESQUEMA IMPLICITO
[F_imp,S_imp,t_imp] = bs_implicito(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

[~,idx_imp] = min(abs(S_imp - S0));
call_imp = F_imp(idx_imp,1);

%% ESQUEMA CRANK-NICOLSON
[F_cn,S_cn,t_cn] = bs_cn(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

[~,idx_cn] = min(abs(S_cn - S0));
call_cn = F_cn(idx_cn,1);

%% ERRORES
err_exp = abs(call_exp - call_exacta);
err_imp = abs(call_imp - call_exacta);
err_cn  = abs(call_cn  - call_exacta);

%% VISUALIZACIÓN DE RESULTADOS
fprintf('---------------------------------------------\n');
fprintf('Comparacion con blsprice (call europea)\n');
fprintf('---------------------------------------------\n');
fprintf('Precio exacto Black-Scholes : %12.6f\n', call_exacta);
fprintf('\n');
fprintf('Esquema explicito          : %12.6f   Error = %.6e\n', call_exp, err_exp);
fprintf('Esquema implicito          : %12.6f   Error = %.6e\n', call_imp, err_imp);
fprintf('Esquema Crank-Nicolson     : %12.6f   Error = %.6e\n', call_cn,  err_cn);
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
% call_bs_curve = zeros(size(S_cn));
% for k = 1:length(S_cn)
%     [call_bs_curve(k), ~] = blsprice(S_cn(k),K,r,T,sigma,q);
% end
% plot(S_cn, call_bs_curve, '--', 'LineWidth', 1.2);
% 
% xline(S0, ':');
% 
% legend('Explícito','Implícito','Crank-Nicolson','Black-Scholes exacta','Location','best');
% xlabel('S');
% ylabel('Precio de la call en t=0');
% title('Comparación de esquemas numéricos para Black-Scholes');
% grid on;
