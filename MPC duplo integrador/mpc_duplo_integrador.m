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
kMax = 500;  % numero maximo de iteracoes
Ts   = 1e-2; % s -> periodo de amostragem
N    = 10;    % horizonte de predicao
xMax = [1.5 1.5]'; % restricao de estado maxima
xMin = [-2 -2]';   % restricao de estado minimo
uMax = 10;         % restricao de entrada maxima
uMin = -10;       % restricao de entrada minima

%% scopes
uLQR = zeros(1, kMax); % entradas lqr
uMPC = zeros(1, kMax); % entradas mpc
u    = zeros(1, kMax); % ulqr+umpc  
x    = zeros(2, kMax); % estados
XMax = vertconcat(xMax, N);
XMin = vertconcat(xMin, N);
UMax = vertconcat(uMax, N);
UMin = vertconcat(uMin, N);
x(:, 1) = x0;          % passando condicao inicial 
options =  optimset('Display','off'); % desabilita logs de quadprog

%% obtendo K do LQR
QLQR = [100 0; 0 10]; % matriz de custo dos estados
RLQR = 1;             % matriz de custo da entrada

% matrizes do sistema
[sysc, sysd] = system_data(Ts);

% calculo do vetor ganhos
[KLQR, S, e] = dlqr(sysd.A, sysd.B, QLQR, RLQR);


%% obtendo matrizes do MPC
QMPC = [1e6 0; 0 5e3]; % matriz de custo dos estados
RMPC = 1;          % matriz de custo da entrada
[A_, Hqp, fqp, Hx, fx, Hu, fu_] = mpc_matrices(sysd.A, sysd.B, QMPC, RMPC, KLQR, N);


%% matrizes de restricao dos estados, controle e terminal
% matrizes que respeitam as restricoes do conjunto Oinf
[Aoinf, boinf, MASObject] = MAS(A_, Ts, KLQR, xMax, xMin, uMax, uMin);

Hab = A_^(N-1)*sysd.B;
for i = N-2:-1:0
    Hab = horzcat(Hab, A_^i*sysd.B);
end

Aoqp = Aoinf*Hab; % restricao terminal
Axqp = [Hx; -Hx]; % restricao de estados
Auqp = [Hu; -Hu]; % restricao de controle

% concatenando todas as restricoes
Aqp = [Axqp; Auqp; Aoqp];

%% definicao da regiao de factibilidade
Af  = [fx; -fx; fu_; -fu_; Aoinf*A_^N];
bxu = [XMax; -XMin; UMax; -UMin; boinf];
ThetaAux = Polyhedron('H', [Aqp Af bxu]);
Theta = ThetaAux.projection(N+1:N+2);

%% simulando acao de controle
for k = 1:kMax
    % calculo acao de controle do LQR
    uLQR(k) = -KLQR*x(:, k); 
    
    % calculo acao de controle MPC
    fqp_ = 2*x(:, k)'*fqp;
    bxqp = [XMax; -XMin] + [-fx; fx]*x(:, k);
    buqp = [UMax; -UMin] + [-fu_; fu_]*x(:, k);
    boqp = boinf - Aoinf*A_^N*x(:, k);
    
    bqp  = [bxqp; buqp; boqp]; % concatena todas as restricoes
    
    UMPC = quadprog(2*Hqp, fqp_, Aqp, bqp, [], [], [], [], [], options);
    uMPC(k) = UMPC(1);
        
    % calculo da acao de controle final
    u(k) = uLQR(k) + uMPC(k);
    
    % aplicando acao de controle
    x(:, k+1) = sysd.A*x(:, k) + sysd.B*u(k); 
    
%     % evoluindo dinamica da planta
%     [t, dummy] = ode45(@(t, x) din_plant(t, x, u(k)), [0 Ts], x(:, k));
%     
%     % atualizando estados
%     x(:, k+1) = dummy(end, :)';
end

%% plotando resultados
% grafico de Theta
figure
Theta.plot('Color',[0,0.7,0.9])

% grafico do MAS
hold on
MASObject.plot()

% plotando estados no grafico do MAS
plot(x(1,:), x(2,:),'y*')

% estados
hold off
figure
plot(linspace(0, kMax*Ts, kMax+1), x)
%stairs(linspace(0, kMax*Ts, kMax+1), x')
ylabel('estados'), xlabel('t[s]'), grid on
legend('x1',  'x2')
title('Duplo integrador')

% entradas
figure
%plot(linspace(0, kMax*Ts, kMax), [u; ulqr; umpc])
stairs(linspace(0, kMax*Ts, kMax), [u; uLQR; uMPC]')
ylabel('entradas'), xlabel('t[s]'), grid on
legend('u = ulqr + umpc', 'ulqr',  'umpc')
title('Duplo integrador')