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
x0   = [0 6*pi/180 0 0]'; % theta, psi, dotTheta, dotPsi
kMax = 500;  % numero maximo de iteracoes da simulacao
Ts   = 4e-3; % s -> periodo de amostragem
N    = 3;    % horizonte de predicao
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
QLQR = [4.2e2 0   0 0
        0   6.5e5 0 0
        0   0   1 0
        0   0   0 1];

% matriz de custo das entradas
RLQR = 1e3*eye(2);

% matrizes do sistema
[sysc, sysd] = system_data(Ts);

% calculo do vetor ganhos
[KLQR, S, e] = dlqr(sysd.A, sysd.B, QLQR, RLQR);


%% obtendo matrizes do MPC
QMPC = [1e2 0     0 0
        0 4.5e6 0 0
        0 0     1 0
        0 0     0 1];
RMPC = 3.5*eye(2);

[Phi, Hqp, fqp, Hx, Px, Hu, Pu] = mpc_matrices(sysd.A, sysd.B, QMPC, RMPC, KLQR, N);

A = sysd.A;
B = sysd.B;

Phi_  = A - B*KLQR;
Px_   = [Phi_; Phi_^2; Phi_^3];
Hx_   = [B        zeros(4, 2) zeros(4, 2)
         Phi_*B   B           zeros(4, 2)
         Phi_^2*B Phi_*B      B          ];
Plqr_ = [-KLQR; zeros(2, 4); zeros(2, 4)];
Hlqr_ = [zeros(2, 4) zeros(2, 4) zeros(2, 4)
         -KLQR       zeros(2, 4) zeros(2, 4)
         zeros(2, 4) -KLQR       zeros(2, 4)];
Pu_ = Plqr_ + Hlqr_*Px_;
Hu_ = Hlqr_*Hx_;
Hu_ = Hu_ + eye(size(Hu_));
diagQ = blkdiag(QMPC, QMPC, QMPC);
diagR = blkdiag(RMPC, RMPC, RMPC);
Hqp_ = Hx_'*diagQ*Hx_ + Hu_'*diagR*Hu_;
fqp_ = Px_'*diagQ*Hx_ + Pu_'*diagR*Hu_;


Fqp_ = 2*x(:, 1)'*fqp_;
[UMPC, fval, exitflag] = quadprog(2*Hqp_, Fqp_);