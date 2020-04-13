function [Phi, Hqp, fqp, Hx, Px, Hu, Pu] = mpc_matrices(A, B, Q, R, K, N)
    % [Hqp, fqp] = mpc_matrices(A, B, Q, R, K, N)
    %
    
    %%
    Phi = A - B*K;
    [m, n] = size(Phi);
    Px = zeros(m*N, n);
    exp = 1;
    for i = 1:m:m*N
        Px(i:i+m-1, 1:n) = Phi^exp;
        exp = exp + 1;
    end
    
    %%
    [m, n] = size(B);
    Hx = zeros(m*N, n*N);
    Hx(1:m, 1:n) = B;
    % preenchendo matriz com Phi^exp * B
    expr = 2;
    for i = m+1:m:m*N
        expc = 1;
        for j = 1:n:n*N
            Hx(i:i+m-1, j:j+n-1) = Phi^(expr-expc)*B;
            expc = expc + 1;
            if (expr-expc) < 0
                break
            end
        end
        expr = expr + 1;
    end
    
    %%
    [m, n] = size(K);
    Plqr = zeros(m*N, n);
    Plqr(1:m, 1:n) = -K;
    Hlqr = zeros(m*N, n*N);
    j = 1;
    for i = m+1:m:m*N
        Hlqr(i:i+m-1, j:j+n-1) = -K;
        j = j + n;
    end
    Pu = Plqr + Hlqr*Px;
    Hu  = Hlqr*Hx + eye(m*N);
    
    %%
    [m, n] = size(Q);
    Q_ = zeros(m*N, n*N);
    j = 1;
    for i = 1:m:m*N
        Q_(i:i+m-1, j:j+n-1) = Q;
        j = j + n;
    end
    
    %%
    [m, n] = size(R);
    R_ = zeros(m*N, n*N);
    j = 1;
    for i = 1:m:m*N
        R_(i:i+m-1, j:j+n-1) = R;
        j = j + n;
    end
    
    %%
    Hqp = Hx'*Q_*Hx + Hu'*R_*Hu;
    fqp = Px'*Q_*Hx + Pu'*R_*Hu;
end