function x_dot = din_plant(t, x, u)
    % Funcao da dinamica nao linear do pendulo de duas rodas.
    % link: https://www.mathworks.com/matlabcentral/fileexchange/19147-nxtway-gs-self-balancing-two-wheeled-robot-controller-design
    %
    % x_dot = din_plant(t, x, u)
    % x_dot -> vetor (4x1) de derivadas dos estados
    %     t -> tempo (nao utilizado na funcao)
    %     x -> vetor (4x1) de estados
    %     u -> vetor (1x2) de entrada
    
    %% parametros da planta
   
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
    
    %% constantes a serem utilizadas nas equacoes

    c1 = (2*m + M)*R^2 + 2*n^2*Jm;
    c2 = M*L*R;
    c3 = 2*n^2*Jm;
    c4 = M*g*L;
    c5 = M*L^2 + Jpsi + 2*n^2*Jm;
    
    h1 = (c2*cos(x(2)) - c3)/c5;
    h2 = (c2*cos(x(2)) - c3)/c1;
    
    %% vetor de estados
    
    x_dot = zeros(4, 1);
    
    %x(1) = theta
    %x(2) = psi
    %x(3) = theta ponto
    %x(4) = psi ponto
    
    x_dot(1) = x(3);
    x_dot(2) = x(4);
    x_dot(3) = (1/(c1 - h1^2*c5))*( -2*(beta*(1 + h1) + fw)*x(3) + 2*beta*(1 + h1)*x(4) ...
                                    + (c2*x(4)^2 - h1*c4)*sin(x(2)) + alpha*(1 + h1)*(u(1) + u(2)) );
    x_dot(4) = (1/(c5 - h2^2*c1))*( 2*(beta + h2*(beta + fw))*x(3) - 2*beta*(1 + h2)*x(4) ...
                                    + (c4 - h2*c2*x(4)^2)*sin(x(2)) - alpha*(1 + h2)*(u(1) + u(2)) );
end