function [Hqp, fqp] = mpc_matrices(Klqr, Q, R, A, B, N)
    
    A_ = A ;
    [m, n] = size(A_);
    Px = zeros(m*N, n);
    exp = 1;
    for i = 1:m:N*m
        Px(i:i+m-1, 1:n) = A_^exp;
        exp = exp + 1;
    end
    
    [m, n] = size(B);
    Hx = zeros(N*m, N*n);
    Hx(1:m, 1:n) = B;
    % preenchendo matriz com A_^exp * B
    expr = 2;
    for i = m+1:m:N*m
        expc = 1;
        for j = 1:n:N*n
            Hx(i:i+m-1, j:j+n-1) = A_^(expr-expc)*B;
            expc = expc + 1;
            if (expr-expc) < 0
                break
            end
        end
        expr = expr + 1;
    end
    
    [m, n] = size(Klqr);
    fu = zeros(m*N, n);
    fu(1:m, 1:n) = -Klqr;
    kappa = zeros(N*m, N*n);
    j = 1;
    for i = m+1:m:N*m
        kappa(i:i+m-1, j:j+n-1) = -Klqr;
        j = j + n;
    end
    
    fu_ = fu + kappa*Px;
    Hu  = kappa*Hx + eye(N);
    
    [m, n] = size(Q);
    Q_ = zeros(N*m, N*n);
    j = 1;
    for i = 1:m:N*m
        Q_(i:i+m-1, j:j+n-1) = Q;
        j = j + n;
    end
    
    [m, n] = size(R);
    R_ = zeros(N*m, N*n);
    j = 1;
    for i = 1:m:N*m
        R_(i:i+m-1, j:j+n-1) = R;
        j = j + n;
    end
    
    Hqp = Hx'*Q_*Hx + Hu'*R_*Hu;
    fqp = Px'*Q_*Hx + fu_'*R_*Hu;
end