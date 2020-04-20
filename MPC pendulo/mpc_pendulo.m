clear
close all
clc

%% implementacao de uma lei de controle preditivo em no modelo linear do pendulo
% Implementacao de uma lei de controle do tipo
%           u(k)    = uLQR(k) + uMPC(k), onde:
%           uLQR(k) = -KLQR*x(k)
%           uMPC(k) = é dado a partir da resolucao de um problema de
%                     otimizacao

%% condicoes iniciais de simulacao
x0   = [0 10*pi/180 0 0]'; % theta, psi, dotTheta, dotPsi
kMax = 1000;  % numero maximo de iteracoes da simulacao
Ts   = 4e-3; % s -> periodo de amostragem
N    = 5;    % horizonte de predicao
xMax = [99999 14*pi/180 99999 99999]'; % restricoes maximas de estados
xMin = -xMax;                      % restricoes minimas de estados
uMax = [9 9]'; % restricoes maximas de controle
uMin = -uMax;      % restricoes minimas de controle

%% scopes
uLQR = zeros(2, kMax); % acoes de controle -> vL, vR
uMPC = zeros(2, kMax);
u    = zeros(2, kMax); % u = uLQR + uMPC
x    = zeros(4, kMax);
XMax = vertconcat(xMax, N);
XMin = vertconcat(xMin, N);
UMax = vertconcat(uMax, N);
UMin = vertconcat(uMin, N);
x(:, 1) = x0; % passando condicao inicial 
options =  optimset('Display','off'); % desabilita logs de quadprog
exitFlags = zeros(1, kMax);           % resultados da funcao quadprog


%% obtendo K do LQR
% matriz de custo dos estados
QLQR = [2.5 0   0 0 % 4.2e2
        0   2.5 0 0 % 6.5e5
        0   0     1 0
        0   0     0 1];

% matriz de custo das entradas
RLQR = 1e-3*eye(2); % 1e3

% matrizes do sistema
[sysc, sysd] = system_data(Ts);

% calculo do vetor ganhos
[KLQR, S, e] = dlqr(sysd.A, sysd.B, QLQR, RLQR);


%% obtendo matrizes do MPC
QMPC = [7.2 0   0 0 % 7.2
        0   4.3 0 0 % 4.3
        0   0     1 0
        0   0     0 1];
RMPC = 4.5*eye(2); % 1.5
[Phi, Hqp, fqp, Hx, Px, Hu, Pu] = mpc_matrices(sysd.A, sysd.B, QMPC, RMPC, KLQR, N);


%% matrizes de restricao dos estados, controle e terminal
% matrizes que respeitam as restricoes do conjunto Oinf
[Aoinf, boinf, MASObject] = MAS(Phi, Ts, KLQR, xMax, xMin, uMax, uMin);

Hab = Phi^(N-1)*sysd.B;
for i = N-2:-1:0
    Hab = horzcat(Hab, Phi^i*sysd.B);
end

Aoqp = Aoinf*Hab; % restricao terminal
Axqp = [Hx; -Hx]; % restricao de estados
Auqp = [Hu; -Hu]; % restricao de controle

% concatenando todas as restricoes
Aqp = [Axqp; Auqp; Aoqp];
%Aqp = [Axqp; Auqp];

%% simulando acoes de controle
for k = 1:kMax
    % calculo acao de controle do LQR
    uLQR(:, k) = -KLQR*x(:, k); 
    
    % calculo acao de controle MPC
    fqp_ = 2*x(:, k)'*fqp;
    bxqp = [XMax; -XMin] + [-Px; Px]*x(:, k);
    buqp = [UMax; -UMin] + [-Pu; Pu]*x(:, k);
    boqp = boinf - Aoinf*Phi^N*x(:, k);
    
    bqp  = [bxqp; buqp; boqp]; % concatena todas as restricoes
    %bqp  = [bxqp; buqp];
    % calcula acoes de controle otimizada
    [UMPC, fval, exitflag] = quadprog(2*Hqp, fqp_, Aqp, bqp, [], [], [], [], [], options);
    if ~isempty(UMPC)
        uMPC(:, k) = UMPC(1:2);
    end
    exitFlags(k) = exitflag;
        
    % calculo da acao de controle final
    u(:, k) = uLQR(:, k) + uMPC(:, k);
    
    % aplicando acao de controle
    x(:, k+1) = sysd.A*x(:, k) + sysd.B*u(:, k); 
    
    fprintf('Iteration: %d\n', k);
end

% estados
figure
plot(linspace(0, kMax*Ts, kMax+1), x)
ylabel('estados'), xlabel('t[s]'), grid on
legend('x1',  'x2', 'x3', 'x4')
title('Pendulo Invertido - Modelo Linear')

% entradas
figure
stairs(linspace(0, kMax*Ts, kMax), [u(1, :); uLQR(1, :); uMPC(1, :)]')
ylabel('entradas'), xlabel('t[s]'), grid on
legend('u = uLQR + uMPC', 'uLQR',  'uMPC')
title('Pendulo invertido linar')

% flags do quadprog
figure
scatter(1:kMax, exitFlags, '.')
ylabel('Quadprog flags'), xlabel('k'), grid on