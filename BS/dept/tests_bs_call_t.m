%% COMPARACION DE ESQUEMAS DE DIFERENCIAS FINITAS PARA BLACK-SCHOLES
% Este script compara los esquemas:
%   - explícito
%   - implícito
%   - Crank-Nicolson
% con una función análoga a blsprice para coeficientes dependientes de t
% en el caso de una call europea vanilla.
% Se ha escrito blsprice_tdep como generalización de la fórmula exacta para
% opciones europeas vanilla cuando r y q dependen de t.

clear; clc; close all;

%% PARAMETROS DEL PROBLEMA
K = 100;      % strike
T = 1.0;      % vencimiento

S0   = 100;   % precio spot donde queremos comparar
Smax = 300;   % truncación del dominio espacial

N = 300;      % número de subintervalos espaciales
M = 4000;     % número de subintervalos temporales

%% COEFICIENTES DEPENDIENTES DEL TIEMPO
% (He puesto ejemplos lineales)
r     = @(t) 0.04 + 0.01*t; 
q     = @(t) 0.01 + 0*t;
sigma = @(t) 0.18 + 0.04*t;

%% DEFINICION DEL PROBLEMA
% Dominio espacial [a,b]
a = 0;
b = Smax;

% Payoff terminal
Phi = @(S) max(S-K,0);

% Condiciones de contorno para call europea
g1 = @(t) 0;
g2 = @(t) Smax*exp(-integral(q,t,T)) - K*exp(-integral(r,t,T)); % Versión análoga obtenida al considerar r y q dep. de t.

%% PRECIO "EXACTO" CON FUNCION ANALOGA A BLSPRICE
[call_exacta, ~] = bsprice_tdep(S0,K,T,r,q,sigma);

%% ESQUEMA EXPLICITO
[F_exp,S_exp,t_exp] = bs_explicito_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

[~,idx_exp] = min(abs(S_exp - S0));
call_exp = F_exp(idx_exp,1);

%% ESQUEMA IMPLICITO
[F_imp,S_imp,t_imp] = bs_implicito_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

[~,idx_imp] = min(abs(S_imp - S0));
call_imp = F_imp(idx_imp,1);

%% ESQUEMA CRANK-NICOLSON
[F_cn,S_cn,t_cn] = bs_cn_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

[~,idx_cn] = min(abs(S_cn - S0));
call_cn = F_cn(idx_cn,1);

%% ERRORES
err_exp = abs(call_exp - call_exacta);
err_imp = abs(call_imp - call_exacta);
err_cn  = abs(call_cn  - call_exacta);

%% VISUALIZACIÓN DE RESULTADOS
fprintf('---------------------------------------------------------\n');
fprintf('Comparacion con bsprice_tdep (call europea, coef. en t)\n');
fprintf('---------------------------------------------------------\n');
fprintf('Precio exacto Black-Scholes : %12.6f\n', call_exacta);
fprintf('\n');
fprintf('Esquema explicito          : %12.6f   Error = %.6e\n', call_exp, err_exp);
fprintf('Esquema implicito          : %12.6f   Error = %.6e\n', call_imp, err_imp);
fprintf('Esquema Crank-Nicolson     : %12.6f   Error = %.6e\n', call_cn,  err_cn);
fprintf('---------------------------------------------------------\n');

%% GRAFICA DE LA SOLUCION EN t=0
% No es especialmente útil pero puede servir para visualizar mejor algunos
% resultados.
% figure;
% plot(S_exp, F_exp(:,1), 'LineWidth', 1.2); hold on;
% plot(S_imp, F_imp(:,1), 'LineWidth', 1.2);
% plot(S_cn , F_cn(:,1) , 'LineWidth', 1.2);
%
% % Fórmula exacta sobre la malla
% call_bs_curve = zeros(size(S_cn));
% for k = 1:length(S_cn)
%     [call_bs_curve(k), ~] = bsprice_tdep(S_cn(k),K,T,r,q,sigma);
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

%% FUNCION LOCAL ANALOGA A BLSPRICE PARA COEFICIENTES DEPENDIENTES DE t
function [Call,Put] = bsprice_tdep(S,K,T,r,q,sigma)
    R = integral(r,0,T);
    Q = integral(q,0,T);
    V = integral(@(u) sigma(u).^2,0,T);

    if V == 0
        Call = max(S*exp(-Q) - K*exp(-R), 0);
        Put  = max(K*exp(-R) - S*exp(-Q), 0);
        return
    end

    d1 = (log(S/K) + R - Q + 0.5*V) / sqrt(V);
    d2 = d1 - sqrt(V);

    Call = S*exp(-Q)*normcdf(d1) - K*exp(-R)*normcdf(d2);
    Put  = K*exp(-R)*normcdf(-d2) - S*exp(-Q)*normcdf(-d1);
end