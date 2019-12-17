clear
close all
clc

%% implementacao de uma lei de controle preditivo em no modelo linear do pendulo
% Implementacao de uma lei de controle do tipo
%           u(k)    = ulqr(k) + umpc(k), onde:
%           ulqr(k) = -Klqr*x(k)
%           umpc(k) = ....


%% condicoes iniciais de simulacao
x0   = [0 5*pi/180 0 0]'; % theta, psi, theta_dot, psi_dot
kmax = 500; % numero maximo de iteracoes da simulacao
Ts   = 4e-3; % s -> periodo de amostragem
N    = 3;    % horizonte de predicao

%% scopes
% linear
ulqr = zeros(2, kmax); % acao de controle -> vl, vr
umpc = zeros(2, kmax); % entradas mpc
u    = zeros(2, kmax); % ulqr+umpc  
x    = zeros(4, kmax); % estados -> theta, psi, theta_dot, psi_dot
x(:, 1) = x0;
options =  optimset('Display','off'); % desabilita logs de quadprog

%% obtendo K do dLQR
% matriz de custo dos estados
Qlqr = [2.4674 0      0 0
        0      2.4674 0 0
        0      0      1 0
        0      0      0 2.4674];

% matriz de custo das entradas
Rlqr = 0.0156*eye(2);

% obtendo matrizes do sistema
[sysc, sysd] = system_data(Ts);

% calculo do ganho  K por meio do dLQR
[Klqr, S, e] = dlqr(sysd.A, sysd.B, Qlqr, Rlqr);

%% obtendo matrizes do MPC
Qmpc = [2.4674 0     0 0
        0     4.5e8 0 0
        0     0     1e3 0
        0     0     0 1e2];
Rmpc = 0.00156*eye(2);
[Hqp, fqp] = mpc_matrices(sysd.A, sysd.B, Qmpc, Rmpc, Klqr, N);

%% simulando acao de controle
for k = 1:kmax
    % calculo acao de controle do LQR
    ulqr(:, k) = -Klqr*x(:, k);
    
    % calculo acao de controle MPC
    fqp_ = 2*x(:, k)'*fqp;
    umpc_aux = quadprog(2*Hqp, fqp_, [], [], [], [], [], [], [], options);
    umpc(:, k) = umpc_aux(1);
    
    % calculo acao de controle final
    u(:, k) = ulqr(:, k) + umpc(:, k);
    
    % aplicando acao de controle
    x(:, k+1) = sysd.A*x(:, k) + sysd.B*u(:, k); 
end


%% plotando resultados
% modelo linear
figure(1)
plot(linspace(0, kmax*Ts, kmax+1), x, 'LineWidth', 2)
ylabel('estados'), xlabel('t[s]'), grid on
legend('theta',  'psi', 'theta ponto', 'psi ponto')
title('Modelo linear')


% entradas
figure(2)
plot(linspace(0, kmax*Ts, kmax), [u; ulqr; umpc])
ylabel('entradas'), xlabel('t[s]'), grid on
legend('u = ulqr + umpc', 'ulqr',  'umpc')
title('Duplo integrador')