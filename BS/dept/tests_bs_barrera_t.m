clear; clc;

% ==========================================================
% Valoracion de opciones barrera con las funciones definidas
% ==========================================================
% Este script reutiliza los esquemas:
%   - bs_explicito_t.m
%   - bs_implicito_t.m
%   - bs_cn_t.m
%
% Permite elegir:
%   - down / up
%   - call / put
%   - knock-out / knock-in
%
% La contraparte se obtiene mediante:
%   V_in + V_out = V_vanilla
% ==========================================================

% ----------------------------------------------------------
% 1. Seleccion del producto
% ----------------------------------------------------------

tipo_barrera = 'down';      % 'down' o 'up'
tipo_opcion  = 'call';    % 'call' o 'put'
tipo_knock   = 'in';      % 'out' o 'in'

metodo = 'cn';            % 'explicito', 'implicito' o 'cn'

% ----------------------------------------------------------
% 2. Parametros financieros
% ----------------------------------------------------------

K  = 100;
T  = 1;
S0 = 100;

r     = @(t) 0.05 + 0*t;
q     = @(t) 0.00 + 0*t;
sigma = @(t) 0.20 + 0*t;

% Barrera
% Para down, S_b debe ser menor que S0.
% Para up, S_b debe ser mayor que S0.

S_b = 80;      % ejemplo down
% S_b = 120;       % ejemplo up

% Frontera artificial para la vanilla y opciones down
S_max = 300;

% Parametros numericos
N = 300;
M = 300;

% ----------------------------------------------------------
% 3. Comprobaciones basicas
% ----------------------------------------------------------

if ~ismember(tipo_barrera, {'down','up'})
    error('tipo_barrera debe ser down o up.');
end

if ~ismember(tipo_opcion, {'call','put'})
    error('tipo_opcion debe ser call o put.');
end

if ~ismember(tipo_knock, {'in','out'})
    error('tipo_knock debe ser in u out.');
end

if ~ismember(metodo, {'explicito','implicito','cn'})
    error('metodo debe ser explicito, implicito o cn.');
end

if strcmp(tipo_barrera,'down') && S0 <= S_b
    error('En una opcion down, la barrera debe estar por debajo de S0.');
end

if strcmp(tipo_barrera,'up') && S0 >= S_b
    error('En una opcion up, la barrera debe estar por encima de S0.');
end

% ----------------------------------------------------------
% 4. Payoff y condiciones de contorno de la vanilla
% ----------------------------------------------------------

R = @(t) integral(r,t,T);
Q = @(t) integral(q,t,T);

switch tipo_opcion

    case 'call'

        Phi = @(S) max(S-K,0);

        g1_vanilla = @(t) 0;

        g2_vanilla = @(t) S_max*exp(-Q(t)) - K*exp(-R(t));

    case 'put'

        Phi = @(S) max(K-S,0);

        g1_vanilla = @(t) K*exp(-R(t));

        g2_vanilla = @(t) 0;
end

% ----------------------------------------------------------
% 5. Valoracion de la opcion europea vanilla
% ----------------------------------------------------------

switch metodo

    case 'explicito'

        [F_vanilla,S_vanilla,t_vanilla] = ...
            bs_explicito_t(0,S_max,T,N,M,r,q,sigma,...
                           Phi,g1_vanilla,g2_vanilla);

    case 'implicito'

        [F_vanilla,S_vanilla,t_vanilla] = ...
            bs_implicito_t(0,S_max,T,N,M,r,q,sigma,...
                           Phi,g1_vanilla,g2_vanilla);

    case 'cn'

        [F_vanilla,S_vanilla,t_vanilla] = ...
            bs_cn_t(0,S_max,T,N,M,r,q,sigma,...
                    Phi,g1_vanilla,g2_vanilla);
end

% Forzamos vectores columna
S_vanilla_col = S_vanilla(:);
V_vanilla_t0  = F_vanilla(:,1);
V_vanilla_t0  = V_vanilla_t0(:);

V_vanilla = interp1(S_vanilla_col,V_vanilla_t0,S0,'linear');

% ----------------------------------------------------------
% 6. Dominio y condiciones de contorno de la knock-out
% ----------------------------------------------------------

switch tipo_barrera

    case 'down'

        % Dominio: [S_b, S_max]
        a_bar = S_b;
        b_bar = S_max;

        % En la barrera: la opcion se extingue
        g1_out = @(t) 0;

        % Frontera superior artificial
        switch tipo_opcion

            case 'call'

                g2_out = @(t) S_max*exp(-Q(t)) - K*exp(-R(t));

            case 'put'

                g2_out = @(t) 0;
        end

    case 'up'

        % Dominio: [0, S_b]
        a_bar = 0;
        b_bar = S_b;

        % Frontera inferior artificial
        switch tipo_opcion

            case 'call'

                g1_out = @(t) 0;

            case 'put'

                g1_out = @(t) K*exp(-R(t));
        end

        % En la barrera: la opcion se extingue
        g2_out = @(t) 0;
end

% ----------------------------------------------------------
% 7. Valoracion de la knock-out
% ----------------------------------------------------------

switch metodo

    case 'explicito'

        [F_out,S_out,t_out] = ...
            bs_explicito_t(a_bar,b_bar,T,N,M,r,q,sigma,...
                           Phi,g1_out,g2_out);

    case 'implicito'

        [F_out,S_out,t_out] = ...
            bs_implicito_t(a_bar,b_bar,T,N,M,r,q,sigma,...
                           Phi,g1_out,g2_out);

    case 'cn'

        [F_out,S_out,t_out] = ...
            bs_cn_t(a_bar,b_bar,T,N,M,r,q,sigma,...
                    Phi,g1_out,g2_out);
end

% Forzamos vectores columna
S_out_col = S_out(:);
V_out_t0  = F_out(:,1);
V_out_t0  = V_out_t0(:);

V_out = interp1(S_out_col,V_out_t0,S0,'linear');

% ----------------------------------------------------------
% 8. Valoracion de la knock-in por paridad
% ----------------------------------------------------------

% Para graficar la knock-in sobre la misma malla que la vanilla,
% interpolamos la knock-out sobre la malla de la vanilla.

V_out_interp_t0 = zeros(size(S_vanilla_col));

for i = 1:length(S_vanilla_col)

    S_actual = S_vanilla_col(i);

    if strcmp(tipo_barrera,'down')

        if S_actual <= S_b
            V_out_interp_t0(i) = 0;
        else
            V_out_interp_t0(i) = interp1(S_out_col,V_out_t0,...
                                         S_actual,'linear',0);
        end

    elseif strcmp(tipo_barrera,'up')

        if S_actual >= S_b
            V_out_interp_t0(i) = 0;
        else
            V_out_interp_t0(i) = interp1(S_out_col,V_out_t0,...
                                         S_actual,'linear',0);
        end

    end
end

% Knock-in por paridad:
% V_in = V_vanilla - V_out
V_in_t0 = V_vanilla_t0 - V_out_interp_t0;

% Evitamos pequenos errores numericos negativos
V_in_t0(V_in_t0 < 0 & abs(V_in_t0) < 1e-10) = 0;

V_in = interp1(S_vanilla_col,V_in_t0,S0,'linear');

% Por seguridad numerica
if V_in < 0 && abs(V_in) < 1e-10
    V_in = 0;
end

% ----------------------------------------------------------
% 9. Seleccion y contraparte
% ----------------------------------------------------------

switch tipo_knock

    case 'out'

        V_seleccionada = V_out;
        V_contraparte  = V_in;

        nombre_seleccionada = [tipo_barrera '-and-out ' tipo_opcion];
        nombre_contraparte  = [tipo_barrera '-and-in ' tipo_opcion];

        S_seleccionada = S_out_col;
        V_seleccionada_t0 = V_out_t0;

        S_contraparte = S_vanilla_col;
        V_contraparte_t0 = V_in_t0;

    case 'in'

        V_seleccionada = V_in;
        V_contraparte  = V_out;

        nombre_seleccionada = [tipo_barrera '-and-in ' tipo_opcion];
        nombre_contraparte  = [tipo_barrera '-and-out ' tipo_opcion];

        S_seleccionada = S_vanilla_col;
        V_seleccionada_t0 = V_in_t0;

        S_contraparte = S_out_col;
        V_contraparte_t0 = V_out_t0;
end

% ----------------------------------------------------------
% 10. Resultados
% ----------------------------------------------------------

fprintf('====================================================\n');
fprintf('Valoracion de opciones barrera\n');
fprintf('Metodo: %s\n', metodo);
fprintf('----------------------------------------------------\n');
fprintf('S0 = %.4f\n', S0);
fprintf('K  = %.4f\n', K);
fprintf('T  = %.4f\n', T);
fprintf('S_b = %.4f\n', S_b);
fprintf('----------------------------------------------------\n');
fprintf('Opcion seleccionada: %s\n', nombre_seleccionada);
fprintf('Precio seleccionada: %.10f\n', V_seleccionada);
fprintf('----------------------------------------------------\n');
fprintf('Contraparte:        %s\n', nombre_contraparte);
fprintf('Precio contraparte: %.10f\n', V_contraparte);
fprintf('----------------------------------------------------\n');
fprintf('Vanilla europea:    %.10f\n', V_vanilla);
fprintf('Knock-out:          %.10f\n', V_out);
fprintf('Knock-in:           %.10f\n', V_in);
fprintf('In + out:           %.10f\n', V_in + V_out);
fprintf('Error paridad:      %.3e\n', abs(V_vanilla - V_in - V_out));
fprintf('====================================================\n');

% ----------------------------------------------------------
% 11. Grafica
% ----------------------------------------------------------

figure;

plot(S_vanilla_col,V_vanilla_t0,'LineWidth',1.5);
hold on;

plot(S_seleccionada,V_seleccionada_t0,'LineWidth',1.5);

plot(S_contraparte,V_contraparte_t0,'--','LineWidth',1.2);

xline(S_b,'--','Barrera');

xlabel('S');
ylabel('Precio');
grid on;

legend('Opcion europea vanilla', ...
       nombre_seleccionada, ...
       nombre_contraparte, ...
       'Barrera', ...
       'Location','Best');

title(['Opcion barrera seleccionada: ' nombre_seleccionada]);
