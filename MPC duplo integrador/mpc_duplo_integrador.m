clear
close all
clc

%% implementacao de uma lei de controle preditivo em um duplo integrador
% Implementacao de uma lei de controle do tipo
%           u(k)    = ulqr(k) + umpc(k), onde:
%           ulqr(k) = -Klqr*x(k)
%           umpc(k) = é dado a partir da resolucao de um problema de
%           otimizacao

%% condicoes iniciais de simulacao
x0   = [1.5 -1]'; 
kmax = 500;  % numero maximo de iteracoes
Ts   = 1e-2; % s -> periodo de amostragem
N    = 3;    % horizonte de predicao
xMax = [1.5 1.5]'; % restricao de estado maxima
xMin = [-2 -2]';   % restricao de estado minimo
uMax = 10;         % restricao de entrada maxima
uMin = -10;       % restricao de entrada minima

%% scopes
ulqr = zeros(1, kmax); % entradas lqr
umpc = zeros(1, kmax); % entradas mpc
u    = zeros(1, kmax); % ulqr+lmpc  
x    = zeros(2, kmax); % estados
Xmax = vertcat(xMax, xMax, xMax);
Xmin = vertcat(xMin, xMin, xMin);
Umax = vertcat(uMax, uMax, uMax);
Umin = vertcat(uMin, uMin, uMin);
x(:, 1) = x0;          % passando condicao inicial 
options =  optimset('Display','off'); % desabilita logs de quadprog

%% obtendo K do lqr
Qlqr = [100 0; 0 10]; % matriz de custo dos estados
Rlqr = 1;          % matriz de custo da entrada

% matrizes do sistema
[sysc, sysd] = system_data(Ts);

% calculo do vetor ganhos
[KLqr, S, e] = dlqr(sysd.A, sysd.B, Qlqr, Rlqr);


%% obtendo matrizes do MPC
Qmpc = [1e6 0; 0 5e3]; % matriz de custo dos estados
Rmpc = 1;          % matriz de custo da entrada
[A_, Hqp, fqp, Hx, fx, Hu, fu_] = mpc_matrices(sysd.A, sysd.B, Qmpc, Rmpc, KLqr, N);


%% matrizes de restricao dos estados e controle
Axqp = vertcat(Hx, -Hx);
bxqp = vertcat(Xmax, -Xmin)*0; % o 0 só serve para dimensionar bxqp

Auqp = vertcat(Hu, -Hu);
buqp = vertcat(Umax, -Umin)*0;% o 0 só serve para dimensionar buqp

%% matrizes de restricao terminal
% matrizes que respeitam as restricoes do conjunto Oinf
[Aoinf, boinf] = MAS(sysd.A, sysd.B, Ts, KLqr,xMax, xMin, uMax, uMin);

Hab = [A_^N*sysd.B A_*sysd.B sysd.B]; % para N = 3

Aoqp = Aoinf*Hab;
boqp = boinf*0; % o 0 só serve para dimensionar boqp


Aqp = vertcat(Axqp, Auqp, Aoqp);

%% simulando acao de controle
for k = 1:kmax
    % calculo acao de controle do LQR
    ulqr(k) = -KLqr*x(:, k); 
    
    % calculo acao de controle MPC
    fqp_ = 2*x(:, k)'*fqp;
    bxqp = vertcat(Xmax, -Xmin) + vertcat(-fx, fx)*x(:, k);
    buqp = vertcat(Umax, -Umin) + vertcat(-fu_, fu_)*x(:, k);
    boqp = boinf - Aoinf*A_^N*x(:, k);
    bqp = vertcat(bxqp, buqp, boqp);
    
    umpc_aux = quadprog(2*Hqp, fqp_, Aqp, bqp, [], [], [], [], [], options);
    if ~isempty(umpc_aux)
        umpc(k) = umpc_aux(1);
    end
    
    % calculo acao de controle final
    u(k) = ulqr(k) + umpc(k);
    
    % aplicando acao de controle
    x(:, k+1) = sysd.A*x(:, k) + sysd.B*u(k); 
    
%     % evoluindo dinamica da planta
%     [t, dummy] = ode45(@(t, x) din_plant(t, x, u(k)), [0 Ts], x(:, k));
%     
%     % atualizando estados
%     x(:, k+1) = dummy(end, :)';
end

%% plotando resultados
% plotando estados no grafico do MAS
hold on
plot(x(1,:), x(2,:),'*')

% estados
hold off
figure
%plot(linspace(0, kmax*Ts, kmax+1), x)
stairs(linspace(0, kmax*Ts, kmax+1), x')
ylabel('estados'), xlabel('t[s]'), grid on
legend('x1',  'x2')
title('Duplo integrador')

% entradas
figure
%plot(linspace(0, kmax*Ts, kmax), [u; ulqr; umpc])
stairs(linspace(0, kmax*Ts, kmax), [u; ulqr; umpc]')
ylabel('entradas'), xlabel('t[s]'), grid on
legend('u = ulqr + umpc', 'ulqr',  'umpc')
title('Duplo integrador')