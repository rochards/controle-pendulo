function [Aoinf, boinf, MASObject] = MAS(Phi, Ts, KLQR, xMax, xMin, uMax, uMin)
    % Funcao que retorna as matrizes do MAS (Maximal Admissible Set) Aoinfo 
    % e boinf, onde Aoinfo*x(k+N|k) <= boinf.
    % Sistema discreto representado no espaco de estados: 
    % x(k+1) = A*x(k)+B*u(k)
    %
    % [Aoinf, boinf] = MAS(A, B, KLqr, xMax, xMin, uMax, uMin)
    %          Aoinf ->
    %          boinf ->
    %              A -> matriz de malha fechada: Phi =  A - B*KLQR
    %             Ts -> [s] - tempo de amostragem
    %           KLQR -> vetor de ganhos do LQR
    %           xMax -> vetor coluna de restricoes maximas de estado
    %           xMin ->    "     "   "       "     minimas "     "
    %           uMax -> vetor coluna de restricoes maximas de entrada
    %           uMin -> "       "    "       "     minimas "     "
    %
    
    %% definindo matrizes que respeitam as restricoes impostas
    % Ax*x(k)<=bx
    [m, n] = size(KLQR);
    Ax = [eye(n); -eye(n); -KLQR; KLQR];
    bx = [xMax; -xMin; uMax; -uMin];
    
    X = Polyhedron('H', [Ax bx]);
    system = LTISystem('A', Phi, 'Ts', Ts);
    MASObject = system.invariantSet('X', X);
    
    %% retornando matrizes do MAS
     
    Aoinf = MASObject.A;
    boinf = MASObject.b;
    
end