[X, Q] = meshgrid(linspace(-2, 2), linspace(-2, 2));

figure
streamslice(X, Q, pHx(X, Q), pHq(X, Q), 'arrow');
hold on;
plot(-1, 0, 'r.', 'markersize', 20);
plot(0, 0, 'b.', 'markersize', 20);
xlabel('x');
ylabel('p');
hold off;

function Hx = pHx(x, p)
    Hx = p + x - x.^3;
end

function Hp = pHq(x, p)
    Hp = p .* (3 * x.^2 - 1);
end