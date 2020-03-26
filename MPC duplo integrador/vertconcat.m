function [X] = vertconcat(x, n)
    % Funcao que realiza em impalhamento vertical de um vetor/matriz n
    % vezes
    %
    % [X] = vertconcat(x, N)
    %   X -> retorno do empilhamento
    %   x -> vetor/matriz a ser empilhado
    %   n -> numero de empilhamento
    %
    X = x;
    for i = 1:n-1
        X = vertcat(X, x);
    end
end