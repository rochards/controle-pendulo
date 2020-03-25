function [Aoinf, boinf] = MAS(A, B, Ts, KLqr, xMax, xMin, uMax, uMin)
    % Funcao que retorna as matrizes do MAS (Maximal Admissible Set) Aoinfo 
    % e boinf, onde Aoinfo*x(k+N|k) <= boinf.
    % Sistema discreto representado no espaco de estados: 
    % x(k+1) = A*x(k)+B*u(k)
    %
    % [Aoinf, boinf] = MAS(A, B, KLqr, xMax, xMin, uMax, uMin)
    %          Aoinf ->
    %          boinf ->
    %              A -> matriz discretizada do sistema 
    %              B ->    "        "       "     "
    %             Ts -> [s] - tempo de amostragem
    %           KLqr -> vetor de ganhos do LQR
    %           xMax -> vetor coluna de restricoes maximas de estado
    %           xMin ->    "     "   "       "     minimas "     "
    %           uMax -> vetor coluna de restricoes maximas de entrada
    %           uMin -> "       "    "       "     minimas "     "
    %
    
    %% definindo matrizes que respeitam as restricoes impostas
    % Ax*x(k)<=bx
    Ax = [eye(2); -eye(2); -KLqr; KLqr];
    bx = [xMax; -xMin; uMax; -uMin];
    
    X = Polyhedron('H', [Ax bx]);
    system = LTISystem('A', A - B*KLqr, 'Ts', Ts);
    invSet = system.invariantSet('X',X);
    
    %% retornando matrizes do MAS
    figure
    invSet.plot()
     
    Aoinf = invSet.A;
    boinf = invSet.b;
    
end