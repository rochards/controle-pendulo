clear
close all
clc

%% implementacao de uma lei de controle preditivo em um duplo integrador
% Implementacao de uma lei de controle do tipo
%           u(k)    = ulqr(k) + umpc(k), onde:
%           ulqr(k) = -Klqr*x(k)
%           umpc(k) = ....

%% condicoes iniciais de simulacao
x0   = [1 1]'; 
kmax = 500;  % numero maximo de iteracoes
Ts   = 1e-2; % s -> periodo de amostragem

%% scopes
ulqr = zeros(1, kmax); % entradas lqr
u    = zeros(1, kmax); % ulqr+lmpc  
x    = zeros(2, kmax); % estados
x(:, 1) = x0;           % passando condicao inicial 

%% obtendo K do lqr
Q = [100 0; 0 10]; % matriz de custo dos estados
R = 1;          % matriz de custo da entrada

% matrizes do sistema
[sysc, sysd] = system_data(Ts);

% calculo do vetor ganhos
[Klqr, S, e] = dlqr(sysd.A, sysd.B, Q, R);

%% simulando acao de controle
for k = 1:kmax
    % calculo acao decontrole do lqr
    ulqr(k) = -Klqr*x(:, k); 
    
    % calculo acao de controle final
    u(k) = ulqr(k);
    
    % evoluindo dinamica da planta
    [t, dummy] = ode45(@(t, x) din_plant(t, x, u(k)), [0 Ts], x(:, k));
    
    % atualizando estados
    x(:, k+1) = dummy(end, :)';
end

%% plotando resultados
plot(linspace(0, kmax*Ts, kmax+1), x, 'LineWidth', 2), hold on
ylabel('estados'), xlabel('t[s]'), grid on
legend('x1',  'x2')
title('Duplo integrador')