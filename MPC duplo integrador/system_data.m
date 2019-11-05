function [sysc, sysd] = system_data(Ts)
    % Funcao que retorna as matrizes de estados do sistema nos tempos 
    % continuo e discreto.
    %
    % [sysc, sysd] = system_data(Ts)
    %        sysc  -> sistema no tempo continuo
    %        sysd  -> sistema no tempo discreto
    %          Ts  -> tempo de discretizacao
    
    %% matrizes da dinamica do sistema
    A = [0 1
         0 0];
    B = [0 
         1];
    C = [1 0];
    D = 0;
    
    %% sistema no espaco de estados no tempo continuo
    sysc = ss(A, B, C, D);
    
    %% sistema no espaco de estados no tempo discreto
    sysd = c2d(sysc, Ts, 'zoh');

end