 function x = findZero(f, x0)
    x = x0;
    m = size(x, 1);
    n = size(f(x), 1);
    criterion = 1e-15;
    d = 1e-15;
    
    while true
        Df = zeros(n, m);
        for i = 1:m
            xp = x; xp(i) = x(i) + d;
            xm = x; xm(i) = x(i) - d;
            Df(:, i) = (f(xp)-f(xm)) / (2*d);
        end

        dx = Df \ f(x);
        if any(isnan(dx))
            error('Met singular point!')
        end
        x = x - dx;
        if norm(dx) < criterion
            break
        end
    end
end

