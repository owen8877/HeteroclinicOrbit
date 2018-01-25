function hess = notOrdinaryHessian(f, x0)
    d = 1e-10;
    n = size(x0, 1);
    m = n / 2;
    if m ~= floor(m)
        error('dimension not even!')
    end
    
    hess = zeros(n, n);
    
    for i = 1:m
        for j = 1:m
            dxi = zeros(n, 1);
            dxi(j) = d;
            dxj = zeros(n, 1);
            dxj(i+m) = d;
            hess(i, j) = ((f(x0+dxi+dxj) - f(x0+dxi-dxj)) + (f(x0-dxi-dxj)  - f(x0-dxi+dxj))) / (4*d*d);
        end
    end
    
    for i = 1:m
        dxi = zeros(n, 1);
        dxi(i+m) = d;
        for j = (m+1):n
            dxj = zeros(n, 1);
            dxj(j) = d;
            hess(i, j) = ((f(x0+dxi+dxj) - f(x0+dxi-dxj)) + (f(x0-dxi-dxj)  - f(x0-dxi+dxj))) / (4*d*d);
        end
    end
    
    for i = (m+1):n
        dxi = zeros(n, 1);
        dxi(i-m) = d;
        for j = 1:m
            dxj = zeros(n, 1);
            dxj(j) = d;
            hess(i, j) = - ((f(x0+dxi+dxj) - f(x0+dxi-dxj)) + (f(x0-dxi-dxj)  - f(x0-dxi+dxj))) / (4*d*d);
        end
    end
    
    for i = (m+1):n
        dxi = zeros(n, 1);
        dxi(i-m) = d;
        for j = (m+1):n
            dxj = zeros(n, 1);
            dxj(j) = d;
            hess(i, j) = - ((f(x0+dxi+dxj) - f(x0+dxi-dxj)) + (f(x0-dxi-dxj)  - f(x0-dxi+dxj))) / (4*d*d);
        end
    end
end