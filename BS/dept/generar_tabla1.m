%% TABLA 1 - Comparacion de esquemas con Black-Scholes exacto
% Muestra una figura con la tabla de resultados.
%
% Compara call y put europeas con coeficientes constantes usando:
%   - Esquema explicito
%   - Esquema implicito
%   - Crank-Nicolson
% frente a blsprice.
%
% El error se calcula en norma infinito sobre toda la malla espacial:
%   ||e||_inf = max_i |f_i^0 - f_BS(S_i,0)|

clear; clc; close all;

%% PARAMETROS DEL PROBLEMA

K = 100;
T = 1.0;

S0    = 100;
S_max = 300;

a = 0;
b = S_max;

N = 300;
M = 4000;   % suficientemente grande para que el explicito sea estable

%% COEFICIENTES CONSTANTES

r0     = 0.05;
q0     = 0.00;
sigma0 = 0.20;

r     = @(t) r0 + 0*t;
q     = @(t) q0 + 0*t;
sigma = @(t) sigma0 + 0*t;

R = @(t) integral(r,t,T);
Q = @(t) integral(q,t,T);

%% PRECIOS EXACTOS EN S0 CON BLSPRICE

[call_exacta, put_exacta] = blsprice(S0,K,r0,T,sigma0,q0);

%% =========================================================
%  CALL EUROPEA
% ==========================================================

Phi_call = @(S) max(S-K,0);

g1_call = @(t) 0;
g2_call = @(t) S_max*exp(-Q(t)) - K*exp(-R(t));

% ---------------------
% Esquema explicito
% ---------------------
tic;
[F_exp_call,S_exp_call,~] = bs_explicito_t(a,b,T,N,M,r,q,sigma,Phi_call,g1_call,g2_call);
call_exp = interp1(S_exp_call,F_exp_call(:,1),S0,'linear');

[call_bs_vec_exp,~] = blsprice(S_exp_call,K,r0,T,sigma0,q0);
call_bs_vec_exp = call_bs_vec_exp(:);

% Correccion por si blsprice devuelve NaN en S=0
call_bs_vec_exp(S_exp_call(:)==0) = 0;

err_call_exp_inf = max(abs(F_exp_call(:,1) - call_bs_vec_exp));
t_call_exp = toc;

% ---------------------
% Esquema implicito
% ---------------------
tic;
[F_imp_call,S_imp_call,~] = bs_implicito_t(a,b,T,N,M,r,q,sigma,Phi_call,g1_call,g2_call);
call_imp = interp1(S_imp_call,F_imp_call(:,1),S0,'linear');

[call_bs_vec_imp,~] = blsprice(S_imp_call,K,r0,T,sigma0,q0);
call_bs_vec_imp = call_bs_vec_imp(:);
call_bs_vec_imp(S_imp_call(:)==0) = 0;

err_call_imp_inf = max(abs(F_imp_call(:,1) - call_bs_vec_imp));
t_call_imp = toc;

% ---------------------
% Crank-Nicolson
% ---------------------
tic;
[F_cn_call,S_cn_call,~] = bs_cn_t(a,b,T,N,M,r,q,sigma,Phi_call,g1_call,g2_call);
call_cn = interp1(S_cn_call,F_cn_call(:,1),S0,'linear');

[call_bs_vec_cn,~] = blsprice(S_cn_call,K,r0,T,sigma0,q0);
call_bs_vec_cn = call_bs_vec_cn(:);
call_bs_vec_cn(S_cn_call(:)==0) = 0;

err_call_cn_inf = max(abs(F_cn_call(:,1) - call_bs_vec_cn));
t_call_cn = toc;

%% =========================================================
%  PUT EUROPEA
% ==========================================================

Phi_put = @(S) max(K-S,0);

g1_put = @(t) K*exp(-R(t));
g2_put = @(t) 0;

% ---------------------
% Esquema explicito
% ---------------------
tic;
[F_exp_put,S_exp_put,~] = bs_explicito_t(a,b,T,N,M,r,q,sigma,Phi_put,g1_put,g2_put);
put_exp = interp1(S_exp_put,F_exp_put(:,1),S0,'linear');

[~,put_bs_vec_exp] = blsprice(S_exp_put,K,r0,T,sigma0,q0);
put_bs_vec_exp = put_bs_vec_exp(:);
put_bs_vec_exp(S_exp_put(:)==0) = K*exp(-r0*T);

err_put_exp_inf = max(abs(F_exp_put(:,1) - put_bs_vec_exp));
t_put_exp = toc;

% ---------------------
% Esquema implicito
% ---------------------
tic;
[F_imp_put,S_imp_put,~] = bs_implicito_t(a,b,T,N,M,r,q,sigma,Phi_put,g1_put,g2_put);
put_imp = interp1(S_imp_put,F_imp_put(:,1),S0,'linear');

[~,put_bs_vec_imp] = blsprice(S_imp_put,K,r0,T,sigma0,q0);
put_bs_vec_imp = put_bs_vec_imp(:);
put_bs_vec_imp(S_imp_put(:)==0) = K*exp(-r0*T);

err_put_imp_inf = max(abs(F_imp_put(:,1) - put_bs_vec_imp));
t_put_imp = toc;

% ---------------------
% Crank-Nicolson
% ---------------------
tic;
[F_cn_put,S_cn_put,~] = bs_cn_t(a,b,T,N,M,r,q,sigma,Phi_put,g1_put,g2_put);
put_cn = interp1(S_cn_put,F_cn_put(:,1),S0,'linear');

[~,put_bs_vec_cn] = blsprice(S_cn_put,K,r0,T,sigma0,q0);
put_bs_vec_cn = put_bs_vec_cn(:);
put_bs_vec_cn(S_cn_put(:)==0) = K*exp(-r0*T);

err_put_cn_inf = max(abs(F_cn_put(:,1) - put_bs_vec_cn));
t_put_cn = toc;

%% TABLA DE RESULTADOS

Opcion = {
    'Call';
    'Call';
    'Call';
    'Put';
    'Put';
    'Put'
};

Metodo = {
    'Explícito';
    'Implícito';
    'Crank-Nicolson';
    'Explícito';
    'Implícito';
    'Crank-Nicolson'
};

Precio_DF = [
    call_exp;
    call_imp;
    call_cn;
    put_exp;
    put_imp;
    put_cn
];

Precio_BS = [
    call_exacta;
    call_exacta;
    call_exacta;
    put_exacta;
    put_exacta;
    put_exacta
];

Error_inf = [
    err_call_exp_inf;
    err_call_imp_inf;
    err_call_cn_inf;
    err_put_exp_inf;
    err_put_imp_inf;
    err_put_cn_inf
];

Tiempo_s = [
    t_call_exp;
    t_call_imp;
    t_call_cn;
    t_put_exp;
    t_put_imp;
    t_put_cn
];

tabla1 = table(Opcion,Metodo,Error_inf,Tiempo_s);

tabla1.Properties.VariableNames = {
    'Opcion', ...
    'Metodo', ...
    'Error_inf', ...
    'Tiempo_s'
};

tabla1.Error_inf = round(tabla1.Error_inf,8);
tabla1.Tiempo_s  = round(tabla1.Tiempo_s,6);

disp(tabla1);

%% MOSTRAR TABLA COMO FIGURA

mostrar_tabla_figura(tabla1);

%% FUNCION AUXILIAR LOCAL

function mostrar_tabla_figura(T)

    fig = figure( ...
        'Color','w', ...
        'Position',[100 100 1150 360], ...
        'Name','Tabla de resultados', ...
        'NumberTitle','off');

    ax = axes(fig);
    axis(ax,'off');
    hold(ax,'on');

    nFilas = height(T);
    nCols  = width(T);

    x0 = 0.04;
    x1 = 0.96;

    yTop        = 0.88;
    yHeaderText = 0.82;
    yMidrule    = 0.76;
    yBottom     = 0.14;

    rowH = (yMidrule - yBottom) / nFilas;

    if nCols == 6
        x = [0.07 0.28 0.49 0.64 0.78 0.91];
    else
        x = linspace(x0,x1,nCols);
    end

    line([x0 x1],[yTop yTop], ...
        'Color','k','LineWidth',1.3);

    line([x0 x1],[yMidrule yMidrule], ...
        'Color','k','LineWidth',0.8);

    line([x0 x1],[yBottom yBottom], ...
        'Color','k','LineWidth',1.3);

    if nFilas == 6
        ySep = yMidrule - 3*rowH;
        line([x0 x1],[ySep ySep], ...
            'Color',[0.45 0.45 0.45], ...
            'LineWidth',0.6);
    end

    nombres = T.Properties.VariableNames;
    nombres_mostrar = nombres;

    for j = 1:numel(nombres_mostrar)
        nombres_mostrar{j} = strrep(nombres_mostrar{j},'_',' ');
    end

    nombres_mostrar(strcmp(nombres,'Opcion'))     = {'Opción'};
    nombres_mostrar(strcmp(nombres,'Metodo'))     = {'Método'};
    nombres_mostrar(strcmp(nombres,'Error_inf'))  = {'||e||_\infty'};
    nombres_mostrar(strcmp(nombres,'Tiempo_s'))   = {'Tiempo (s)'};

    for j = 1:nCols
        text(x(j),yHeaderText,nombres_mostrar{j}, ...
            'HorizontalAlignment','center', ...
            'FontSize',12.5, ...
            'FontWeight','bold', ...
            'Interpreter','tex');
    end

    datos = table2cell(T);

    for i = 1:nFilas

        y = yMidrule - (i-0.5)*rowH;

        for j = 1:nCols

            valor = datos{i,j};
            nombre_col = nombres{j};

            if isstring(valor)
                txt = char(valor);

            elseif ischar(valor)
                txt = valor;

            elseif isnumeric(valor)

                if strcmp(nombre_col,'Error_inf')
                    if abs(valor) >= 1e5 || abs(valor) < 1e-4
                        txt = sprintf('%.3e',valor);
                    else
                        txt = sprintf('%.6f',valor);
                    end

                elseif strcmp(nombre_col,'Tiempo_s')
                    txt = sprintf('%.4f',valor);

                else
                    txt = sprintf('%.6g',valor);
                end

            else
                txt = char(string(valor));
            end

            text(x(j),y,txt, ...
                'HorizontalAlignment','center', ...
                'FontSize',11.5, ...
                'Interpreter','none');
        end
    end

    xlim([0 1]);
    ylim([0 1]);

end
