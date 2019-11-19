clear
close all
clc

Ts = 1e-2; % tempo de amostragem

%% matrizes do mpc


%% matrizes do sistema
A = [0 1; 0 0];
B = [0; 1];
C = [1 0];
D = 0;

%% sistema no espaco de estados no tempo continuo
sysc = ss(A, B, C, D);

%% sistema no espaco de estados no tempo discreto
sysd = c2d(sysc, Ts, 'zoh');

%% calculo vetor Klqr

Qlqr = [100 0; 0 10]; % matriz de custo dos estados
Rlqr = 1;             % matriz de custo das entradas

% calculo do vetor ganhos
[Klqr, S, e] = dlqr(sysd.A, sysd.B, Qlqr, Rlqr);

%% montando matrizes do mpc para N = 3

Qmpc = [1e6 0; 0 5e3]; % matriz de custo dos estados
Rmpc = 1;          % matriz de custo da entrada

A_ = sysd.A - sysd.B*Klqr;

fx = vertcat(A_, A_^2, A_^3)

Hx = vertcat(horzcat(sysd.B, [0; 0], [0; 0]),...
             horzcat(A_*sysd.B, sysd.B, [0; 0]),...
             horzcat(A_^2*sysd.B, A_*sysd.B, sysd.B))

fu = vertcat(-Klqr, [0 0], [0 0])

kappa = vertcat(horzcat([0 0], [0 0], [0 0]),...
                horzcat(-Klqr, [0 0], [0 0]),...
                horzcat([0 0], -Klqr, [0 0]))
            

fu_ = fu + kappa*fx

Hu = kappa*Hx + eye(3)

Q_ = blkdiag(Qmpc, Qmpc, Qmpc)
R_ = blkdiag(Rmpc, Rmpc, Rmpc)


Hqp = Hx'*Q_*Hx + Hu'*R_*Hu
fqp = fx'*Q_*Hx + fu_'*R_*Hu 