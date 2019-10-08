clear
close
clc

% implementacao pendulo invertido

function x_dot = din_plant(x, vl, vr)
    % == parametros da planta ==
    g = 9.81;  % m/s� -> aceleracao da gravidade
    m = 0.03;  % kg   -> massa das rodas
    R = 0.04;  % m    -> raio das rodas
    M = 0.6;   % kg   -> massa do corpo (p�ndulo)
    L = 0.072; % m    -> dist�ncia da roda ao centro de massa
    Jpsi = 0.001;  % kgm�     -> body pitch inertia moment
    Jm   = 10^-5;  % kgm�     -> momento de inercia do motor DC
    Rm   = 6.69;   % ohm      -> resistencia do motor DC
    Kb   = 0.468;  % V s/rad  -> DC motor back EMF constant
    Kt   = 0.317;  % Nm/A     -> constante do torque do motor DC
    fm   = 0.0022; % Nm rad/s -> coeficiente de atrito do motor DC
    Vb   = 8.00;   % V        -> tensao de alimentacao
    Vo   = 0.625;  % V        -> tensao de offset
    mi   = 1.089;  %          -> fator de ganho da alimentacao
    Gu   = 10^-2;  %          -> fator de ganho do PWM
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
    
    f1 = (c2*cos(x(2)) - c3)/c5;
    f2 = (c2*cos(x(2)) - c3)/c1;
    
    % == iniciando o vetor de estados == 
    x_dot    = zeros(4, 1);
    x_dot(1) = x(3);
    x_dot(2) = x(4);
    x_dot(3) = (1/(c1 - f1*(c2*cos(x(2)) - c3)))*( -2*(beta*(1 + f1) + fw)*x(3) + 2*beta*(1 + f1)*x(4) ...
                                                  + (c2*x(4)^2 - f1)*sin(x(2)) + alpha*(1 + f1)*(vl + vr) );
    x_dot(4) = (1/(c5 - f2*(c2*cos(x(2)) - c3)))*( 2*(beta + f2*(beta + fw))*x(3) - 2*beta*(1 + f2)*x(4) ...
                                                  + (c4 - f2*c2*x(4)^2)*sin(x(2)) - alpha*(1 + f2)*(vl + vr) );
    
end