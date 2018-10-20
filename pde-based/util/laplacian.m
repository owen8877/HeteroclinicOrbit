function l = laplacian(x, h)
    [m, n] = size(x);
    l = zeros(m, n);
    l(:, 2:end-1) = (x(:, 1:end-2) + x(:, 3:end) - 2*x(:, 2:end-1)) / h^2;
    l(:, 1) = l(:, 2);
    l(:, end) = l(:, end-1);
end