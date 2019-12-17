clear
close all
clc

Ts = 4e-3; % tempo de amostragem

%% matrizes do mpc


%% matrizes do sistema

[sysc, sysd] = system_data(Ts);


%% calculo vetor Klqr

% matriz de custo dos estados
Qlqr = [2.4674 0      0 0
            0      2.4674 0 0
            0      0      1 0
            0      0      0 2.4674];

% matriz de custo das entradas
Rlqr = 0.0156*eye(2);

% calculo do vetor ganhos
[Klqr, S, e] = dlqr(sysd.A, sysd.B, Qlqr, Rlqr);

%% montando matrizes do mpc para N = 3

Qmpc = 1e6*[2.4674 0      0 0
        0      2.4674 0 0
        0      0      1 0
        0      0      0 2.4674];
Rmpc = eye(2);

A_ = sysd.A - sysd.B*Klqr;

fx = vertcat(A_, A_^2, A_^3)

Hx = vertcat(horzcat(sysd.B, zeros(4,2), zeros(4,2)),...
             horzcat(A_*sysd.B, sysd.B, zeros(4,2)),...
             horzcat(A_^2*sysd.B, A_*sysd.B, sysd.B))

fu = vertcat(-Klqr, zeros(2, 4), zeros(2, 4))

kappa = vertcat(horzcat(zeros(2, 4), zeros(2, 4), zeros(2, 4)),...
                horzcat(-Klqr, zeros(2, 4), zeros(2, 4)),...
                horzcat(zeros(2, 4), -Klqr, zeros(2, 4)))
            

fu_ = fu + kappa*fx

Hu = kappa*Hx + eye(6)

Q_ = blkdiag(Qmpc, Qmpc, Qmpc)
R_ = blkdiag(Rmpc, Rmpc, Rmpc)


Hqp = Hx'*Q_*Hx + Hu'*R_*Hu
fqp = fx'*Q_*Hx + fu_'*R_*Hu 