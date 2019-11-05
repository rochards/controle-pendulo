function dot_x = din_plant(t, x, u)
    % Funcao que representa a dinamica de um duplo integrador
    %
    % dot_x = din_plant(t, x, u)
    % dot_x -> vetor (2x1) de derivadas dos estados
    %     t -> tempo (nao utilizado na funcao)
    %     x -> vetor (2x1) de estados
    %     u -> vetor (1x1) de entrada
    
    %% vetor de estados
    dot_x = zeros(2, 1);
    dot_x(1) = x(2);
    dot_x(2) = u;
end