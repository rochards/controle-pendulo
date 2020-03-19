% definição de região de factibilidade
clear
close all
clc

%% condicoes iniciais
Ts = 1e-2; % [s] -> periodo de amostragem
xMax = [1.5 1.5]'; % restricao maxima de estado
xMin = [-2 -2]';   % restricao minima de estado
uMax = 10;   % restricao maxima de entrada
uMin = -10; % restricao minima de entrada

%% obtendo K do lqr
Qlqr = [100 0; 0 10]; % matriz de custo dos estados
Rlqr = 1;          % matriz de custo da entrada

% matrizes do sistema
[sysc, sysd] = system_data(Ts);

% calculo do vetor ganhos
[KLqr, S, e] = dlqr(sysd.A, sysd.B, Qlqr, Rlqr);

%% utilizando o MPT
Ax = [eye(2); -eye(2); -KLqr; KLqr];
bx = [xMax; -xMin; uMax; -uMin];
PolyX = Polyhedron('H', [Ax bx] );
system = LTISystem('A', sysd.A - sysd.B*KLqr);
% system.x.max = x_max;
% system.x.min = x_min;
% system.u.max = u_max;
% system.u.min = u_min;
invSet = system.invariantSet('X',PolyX);
invSet.plot()