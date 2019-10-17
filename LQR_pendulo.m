clear
close
clc

% == implementacao pendulo invertido ==

% == definindo variaveis iniciais ==
x0   = [0; pi/180; 0; 0];   % condicao incial -> theta; psi; theta_dot; psi_dot
imax = 1000;             % numero maximo de iteracoes da simulacao
Ts   = 4e-3;           % s -> periodo de amostragem
y    = zeros(imax, 3); % saída -> theta; theta_dot; psi_dot
u    = zeros(2, imax); % acao de controle -> vl, vr

psi = zeros(imax, 1);

% == calculo da matriz de ganhos ==  
K    = gain_matrix(Ts);

for i = 2:imax
    
    % calculando acao de controle
    u(:,i) = -K*x0;
    
    % == evoluindo a dinamica da planta ==
    [t, x_states] = ode45(@(t, x) din_plant(t, x, u(:,i)), [0 Ts], x0);
    
    % == atualizando estados da planta ==
    x0(1:4) = x_states(end, 1:4);
    
    % == atualizando saidas ==
    y(i+1, 1) = x_states(end, 1); % theta
    psi(i+1)  = x_states(end, 2)*180/pi; 
    y(i+1, 2) = x_states(end, 3); % theta_dot
    y(i+1, 3) = x_states(end, 4); % psi_dot
end

% == plotando os resultados para psi
plot(linspace(0, imax*Ts, imax+1), psi(1:end), 'LineWidth', 2), hold on
ylabel('psi[º]'), xlabel('t[s]'), grid on

function K = gain_matrix(Ts)
    % == parametros da planta ==
    g = 9.81;  % m/s² -> aceleracao da gravidade
    m = 0.03;  % kg   -> massa das rodas
    R = 0.04;  % m    -> raio das rodas
    M = 0.6;   % kg   -> massa do corpo (pêndulo)
    L = 0.072; % m    -> distância da roda ao centro de massa
    Jpsi = 0.001;  % kgm²     -> body pitch inertia moment
    Jm   = 1e-5;   % kgm²     -> momento de inercia do motor DC
    Rm   = 6.69;   % ohm      -> resistencia do motor DC
    Kb   = 0.468;  % V s/rad  -> DC motor back EMF constant
    Kt   = 0.317;  % Nm/A     -> constante do torque do motor DC
    fm   = 0.0022; % Nm rad/s -> coeficiente de atrito do motor DC
    Vb   = 8.00;   % V        -> tensao de alimentacao
    Vo   = 0.625;  % V        -> tensao de offset
    mi   = 1.089;  %          -> fator de ganho da alimentacao
    Gu   = 1e-2;   %          -> fator de ganho do PWM
    n    = 1;      %          -> gearbox ratio
    fw   = 0;      %          -> coeficiente de atrito entre a roda e o chao.

    alpha = n*Kt/Rm;
    beta  = n*Kt*Kb/Rm + fm;
    
    % == definicao das constantes para serem utilizadas nas equacoes ==
    c1 = (2*m + M)*R^2 + 2*n^2*Jm;
    c2 = M*L*R;
    c3 = 2*n^2*Jm;
    c4 = M*g*L;
    c5 = M*L^2 + Jpsi + 2*n^2*Jm;
    
    % == definindo derivadas parciais das matrizes ==
    df1dx1 = 0; df1dx2 = 0; df1dx3 = 1; df1dx4 = 0;
    df2dx1 = 0; df2dx2 = 0; df2dx3 = 0; df2dx4 = 1;
    df3dx1 = 0;
    df3dx2 = -(c2 - c3)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df3dx3 = -2*(beta*(c2 - c3 + c5) + c5*fw)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df3dx4 = 2*beta*(c2 - c3 + c5)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df4dx1 = 0;
    df4dx2 = (c1*c4)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df4dx3 = 2*(beta*(c1 + c2 - c3) + fw*(c2 - c3))/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df4dx4 = -2*beta*(c1 + c2 - c3)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df1du1 = 0; df1du2 = 0;
    df2du1 = 0; df2du2 = 0;
    df3du1 = alpha*(c2 - c3 + c5)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df3du2 = df3du1;
    df4du1 = - alpha*(c1 + c2 - c3)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df4du2 = df4du1;
    
    % == definindo matrizes da dinamica do sistema ==
    A = [df1dx1 df1dx2 df1dx3 df1dx4
         df2dx1 df2dx2 df2dx3 df2dx4
         df3dx1 df3dx2 df3dx3 df3dx4
         df4dx1 df4dx2 df4dx3 df4dx4];
    
    B = [df1du1 df1du2
         df2du1 df2du2
         df3du1 df3du2
         df4du1 df4du2];
  
    C = [1 0 0 0
         0 0 1 0
         0 0 0 1];
     
    % == criando o sistema no espaço de estados no tempo discreto ==
    sys = ss(A, B, C, 0, Ts);
    
    % == implementando o dLQR ==
    Q = [1 0 0 0
         0 6e5 0 0
         0 0 1 0
         0 0 0 1];
    R = eye(2);
    [K, S, e] = dlqr(sys.A, sys.B, Q, R, 0);
end

function x_dot = din_plant(t, x, u)
    % == parametros da planta ==
    g = 9.81;  % m/s² -> aceleracao da gravidade
    m = 0.03;  % kg   -> massa das rodas
    R = 0.04;  % m    -> raio das rodas
    M = 0.6;   % kg   -> massa do corpo (pêndulo)
    L = 0.072; % m    -> distância da roda ao centro de massa
    Jpsi = 0.001;  % kgm²     -> body pitch inertia moment
    Jm   = 1e-5;   % kgm²     -> momento de inercia do motor DC
    Rm   = 6.69;   % ohm      -> resistencia do motor DC
    Kb   = 0.468;  % V s/rad  -> DC motor back EMF constant
    Kt   = 0.317;  % Nm/A     -> constante do torque do motor DC
    fm   = 0.0022; % Nm rad/s -> coeficiente de atrito do motor DC
    Vb   = 8.00;   % V        -> tensao de alimentacao
    Vo   = 0.625;  % V        -> tensao de offset
    mi   = 1.089;  %          -> fator de ganho da alimentacao
    Gu   = 1e-2;   %          -> fator de ganho do PWM
    n    = 1;      %          -> gearbox ratio
    fw   = 0;      %          -> coeficiente de atrito entre a roda e o chao.

    alpha = n*Kt/Rm;
    beta  = n*Kt*Kb/Rm + fm;
    
    % == definicao das constantes para serem utilizadas nas equacoes ==
    c1 = (2*m + M)*R^2 + 2*n^2*Jm;
    c2 = M*L*R;
    c3 = 2*n^2*Jm;
    c4 = M*g*L;
    c5 = M*L^2 + Jpsi + 2*n^2*Jm;
    
    h1 = (c2*cos(x(2)) - c3)/c5;
    h2 = (c2*cos(x(2)) - c3)/c1;
    
    % == tensoes aplicadas aos motores ==
    %vl = Gu*(mi*Vb - Vo)*u; % motor esquerdo
    %vr = vl;                % motor direito
    u(1) = Gu*(mi*Vb - Vo)*u(1);
    u(2) = Gu*(mi*Vb - Vo)*u(2);
    
    % == iniciando o vetor de estados == 
    x_dot    = zeros(4, 1);
    x_dot(1) = x(3);
    x_dot(2) = x(4);
    x_dot(3) = (1/(c1 - h1*(c2*cos(x(2)) - c3)))*( -2*(beta*(1 + h1) + fw)*x(3) + 2*beta*(1 + h1)*x(4) ...
                                                  + (c2*x(4)^2 - h1)*sin(x(2)) + alpha*(1 + h1)*(u(1) + u(2)) );
    x_dot(4) = (1/(c5 - h2*(c2*cos(x(2)) - c3)))*( 2*(beta + h2*(beta + fw))*x(3) - 2*beta*(1 + h2)*x(4) ...
                                                  + (c4 - h2*c2*x(4)^2)*sin(x(2)) - alpha*(1 + h2)*(u(1) + u(2)) );
    
end