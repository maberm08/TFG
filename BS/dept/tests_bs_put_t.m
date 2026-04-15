%% COMPARACION DE ESQUEMAS DE DIFERENCIAS FINITAS PARA BLACK-SCHOLES
% Este script compara los esquemas:
%   - explícito
%   - implícito
%   - Crank-Nicolson
% con una función análoga a blsprice para coeficientes dependientes de t
% en el caso de una put europea vanilla.
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
Phi = @(S) max(K-S,0);

% Condiciones de contorno para put europea
g1 = @(t) K*exp(-integral(r,t,T)); % Versión análoga obtenida al considerar r dep. de t.
g2 = @(t) 0;

%% PRECIO "EXACTO" CON FUNCION ANALOGA A BLSPRICE
[~, put_exacta] = bsprice_tdep(S0,K,T,r,q,sigma);

%% ESQUEMA EXPLICITO
[F_exp,S_exp,t_exp] = bs_explicito_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

[~,idx_exp] = min(abs(S_exp - S0));
put_exp = F_exp(idx_exp,1);

%% ESQUEMA IMPLICITO
[F_imp,S_imp,t_imp] = bs_implicito_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

[~,idx_imp] = min(abs(S_imp - S0));
put_imp = F_imp(idx_imp,1);

%% ESQUEMA CRANK-NICOLSON
[F_cn,S_cn,t_cn] = bs_cn_t(a,b,T,N,M,r,q,sigma,Phi,g1,g2);

[~,idx_cn] = min(abs(S_cn - S0));
put_cn = F_cn(idx_cn,1);

%% ERRORES
err_exp = abs(put_exp - put_exacta);
err_imp = abs(put_imp - put_exacta);
err_cn  = abs(put_cn  - put_exacta);

%% VISUALIZACIÓN DE RESULTADOS
fprintf('---------------------------------------------------------\n');
fprintf('Comparacion con bsprice_tdep (put europea, coef. en t)\n');
fprintf('---------------------------------------------------------\n');
fprintf('Precio exacto Black-Scholes : %12.6f\n', put_exacta);
fprintf('\n');
fprintf('Esquema explicito          : %12.6f   Error = %.6e\n', put_exp, err_exp);
fprintf('Esquema implicito          : %12.6f   Error = %.6e\n', put_imp, err_imp);
fprintf('Esquema Crank-Nicolson     : %12.6f   Error = %.6e\n', put_cn,  err_cn);
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
% put_bs_curve = zeros(size(S_cn));
% for k = 1:length(S_cn)
%     [~, put_bs_curve(k)] = bsprice_tdep(S_cn(k),K,T,r,q,sigma);
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