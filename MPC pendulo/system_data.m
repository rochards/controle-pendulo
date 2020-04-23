function [sysc, sysd] = system_data(Ts)
    % Funcao que retorna as matrizes de estados do sistema nos tempos 
    % continuo  e discreto a partir do modelo linearizado do sistema.
    %
    % [sysc, sysd] = system_data(Ts)
    %        sysc  -> sistema no tempo continuo
    %        sysd  -> sistema no tempo discreto
    %          Ts  -> tempo de discretizacao

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
    
    %% derivadas parciais das matrizes
    
    df1dx1 = 0;
    df1dx2 = 0;
    df1dx3 = 1;
    df1dx4 = 0;
    
    df2dx1 = 0;
    df2dx2 = 0;
    df2dx3 = 0;
    df2dx4 = 1;
    
    df3dx1 = 0;
    df3dx2 = -c4*(c2 - c3)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df3dx3 = -2*(beta*(c2 - c3 + c5) + c5*fw)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df3dx4 = 2*beta*(c2 - c3 + c5)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    
    df4dx1 = 0;
    df4dx2 = c1*c4/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df4dx3 = 2*(beta*(c1 + c2 - c3) + fw*(c2 - c3))/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df4dx4 = -2*beta*(c1 + c2 - c3)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    
    df1du1 = 0;
    df1du2 = 0;
    
    df2du1 = 0;
    df2du2 = 0;
    
    df3du1 = alpha*(c2 - c3 + c5)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);%ok
    df3du2 = df3du1;
    
    df4du1 = -alpha*(c1 + c2 - c3)/(-c2^2 + 2*c2*c3 - c3^2 + c1*c5);
    df4du2 = df4du1;
    
    %% matrizes da dinamica do sistema
    
    A = [df1dx1 df1dx2 df1dx3 df1dx4
         df2dx1 df2dx2 df2dx3 df2dx4
         df3dx1 df3dx2 df3dx3 df3dx4
         df4dx1 df4dx2 df4dx3 df4dx4];
    
    B = [df1du1 df1du2
         df2du1 df2du2
         df3du1 df3du2
         df4du1 df4du2];
    
   %C = [1 0 0 0; 0 1 0 0];
   C = eye(4); 
   
   D = 0;
   
   %% sistema no espaco de estados no tempo continuo
   
   sysc = ss(A, B, C, D);
   
   %% sistema no espaco de estados no tempo discreto
   
   sysd = c2d(sysc, Ts, 'zoh');
end