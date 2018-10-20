function d = derivative(v, h)
    [m, n] = size(v);
    d = zeros(m, n);
    d(:, 1:end-1) = (v(:, 2:end)-v(:, 1:end-1))/h;
    d(:, end) = d(:, end-1);
end