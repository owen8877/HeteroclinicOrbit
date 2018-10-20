function K = newtonSubSolver(x, phi)
    K = phi(x);
    for itr = 1:20
        nK = phi(K);
        if norm(nK-K) < 1e-12
            K = nK;
            break
        end
        K = nK;
    end
    if itr == 20
        fprintf('Maybe not accurate.\n')
    end
    return

    count = 0;
    n = numel(x);
    Jf = zeros(n);
    dh = 1e-8;
    for i = 1:n
        e = zeros(n, 1);
        e(i) = dh;
        Jf(:, i) = (phi(x+e) - phi(x)) / dh;
    end
    J = eye(n) - Jf;
    K = x;
    
    while true
        f = K - phi(K);
        
        if norm(f) < 1e-15
            break
        end
        
        if count ~= 0
            df = f - oldf;
            dK = K - oldK;
            J = J + (df - J*dK) / (norm(dK)^2) * dK';
        end
        
        oldK = K;
        oldf = f;
        K = K - J \ f;
    end
end