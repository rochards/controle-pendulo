clear
close all
clc

%% implementacao pendulo invertido

% Implementacao de uma lei de controle u(k) = -K*x(k), onde o vetor de 
% ganhos K eh obtido pelo dLQR.
% As simulacoes foram realizadas para os modelos linear e nao linear do 
% pendulo.
% A linearizacao foi feita por serie de Taylor.


%% condicoes iniciais de simulacao

x0   = [0 30*pi/180 0 0]'; % theta, psi, theta_dot, psi_dot
kmax = 2000; % numero maximo de iteracoes da simulacao
Ts   = 4e-3; % s -> periodo de amostragem
options2 = odeset('Reltol', 1e-7,'AbsTol', 1e-7); %ode45

%% scopes

% nao linear
unl = zeros(2, kmax); % acao de controle -> vl, vr
xnl = zeros(4, kmax); % estados -> theta, psi, theta_dot, psi_dot
xnl(:, 1) = x0;

% linear
ul = zeros(2, kmax); % acao de controle -> vl, vr
xl = zeros(4, kmax); % estados -> theta, psi, theta_dot, psi_dot
xl(:, 1) = x0;

%% obtencao do vetor de ganhos K

% matriz de custo dos estados
Q = [4e2 0      0 0
     0      6e5 0 0
     0      0      1 0
     0      0      0 1];

% matriz de custo das entradas
R = 1e3*eye(2);

% obtendo matrizes do sistema
[sysc, sysd] = system_data(Ts);

% calculo do ganho  K por meio do dLQR
[K, S, e] = dlqr(sysd.A, sysd.B, Q, R);

%% simulando acao de controle

for k = 1:kmax
    
    % == modelo nao linear ==
    
    unl(:, k) = -K*xnl(:, k); % calculando acao de controle
    [t, x] = ode45(@(t, x) din_plant(t, x, unl(:, k)), [0 Ts], xnl(:,k)); % evoluindo a dinamica da planta
    xnl(:, k+1) = x(end, :)'; % atualizando estados da planta
    
    
    % == modelo linear == 
    
    ul(:, k) = -K*xl(:, k); % calculando acao de controle
    xl(:, k+1) = sysd.A*xl(:, k) + sysd.B*ul(:, k); 
    
end

%% plotando resultados

% modelo nao linear
figure(1)
subplot(2, 1, 1)
plot(linspace(0, kmax*Ts, kmax+1), xnl, 'LineWidth', 2), hold on
ylabel('estados'), xlabel('t[s]'), grid on
legend('theta',  'psi', 'theta ponto', 'psi ponto')
title('Modelo nao linear')

% modelo linear
subplot(2, 1, 2)
plot(linspace(0, kmax*Ts, kmax+1), xl, 'LineWidth', 2), hold on
ylabel('estados'), xlabel('t[s]'), grid on
legend('theta',  'psi', 'theta ponto', 'psi ponto')
title('Modelo linear')

% valores de psi para ambos os modelos
figure(2)
plot(linspace(0, kmax*Ts, kmax+1), [xnl(2, :); xl(2, :)], 'LineWidth', 2), hold on
ylabel('psi [rad]'), xlabel('t[s]'), grid on
legend('nao linear',  'linear')
title('psi')

% entradas
figure(3)
plot(linspace(0, kmax*Ts, kmax), ul(1,:))
title('modelo linear')
ylabel('entrada'), xlabel('t[s]'), grid on