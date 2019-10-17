clear
close
clc

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

Jw    = m*R^2/2; % Wheel inertia moment
alpha = n*Kt/Rm;
beta  = n*Kt*Kb/Rm + fm;



E = [((2*m + M)*R^2 + 2*Jw + 2*n^2*Jm) (M*L*R - 2*n^2*Jm)
     (M*L*R - 2*n^2*Jm)                (M*L^2 + Jpsi + 2*n^2*Jm)];
F = 2*[(beta + fw) -beta
        -beta       beta];
G = [0, 0; 0 -M*g*L];
H = alpha*[1, 1; -1, -1];


A = zeros(4, 4);
A(1, 3) = 1;
A(2, 4) = 1;
A(3, 2) = -g*M*L*E(1,2)/det(E);
A(3, 3) = -2*((beta + fw)*E(2,2) + beta*E(1,2))/det(E);
A(3, 4) = 2*beta*(E(2,2) + E(1,2))/det(E);
A(4, 2) = g*M*L*E(1,1)/det(E);
A(4, 3) = 2*((beta + fw)*E(1,2) + beta*E(1,1))/det(E);
A(4, 4) = -2*beta*(E(1,1) + E(1,2))/det(E);

B = zeros(2, 2);
B(3, 1) = alpha*(E(2,2) + E(1,2))/det(E);
B(3, 2) = B(3, 1);
B(4, 1) = -alpha*(E(1,1) + E(1,2))/det(E);
B(4, 2) = B(4, 1);

C = [1 0 0 0
     0 1 0 0
     0 0 1 0
     0 0 0 1];

 
Ts = 4e-3;
sys = ss(A, B, C, 0, Ts);
    
Q = [1 0 0 0
     0 1 0 0
     0 0 1 0
     0 0 0 1];
R = [1 0
     0 1];
[K, S, e] = dlqr(sys.A, sys.B, Q, R, 0);

x0 = [0; 5*pi/180; 0; 0];

sys = ss(A - B*K, B, C, 0, Ts);
t = 0:4e-3:4;
[y, t, x] = initial(sys, x0, t);

plot(t, y(:,2))