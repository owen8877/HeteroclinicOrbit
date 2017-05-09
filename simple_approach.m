[X, Q] = meshgrid(linspace(-2, 1), linspace(-3, 1));

figure
streamslice(X, Q, pHx(X, Q), pHq(X, Q), 'arrow');
hold on;
plot(0, 0, 'r.', 'markersize', 20);
plot(-1, -2, 'b.', 'markersize', 20);
xlabel('x');
ylabel('q');
hold off;

tMax = 5;
tPositiveSpan = [0 tMax];
tNegativeSpan = [tMax 0];
start = 1e-7;

[t1p, q1] = ode45(@odeFunctionRound1Positive, tPositiveSpan, -2 + start);
solutionRound1Positive = ode45(@odeFunctionRound1Positive, tPositiveSpan, -2 + start);
q1_overflow_index = find(q1 > 0);
q1 = q1([1 q1_overflow_index(1)]);
t1p = t1p([1 q1_overflow_index(1)]);
figure
plot(t1p, q1, 'r-o')

[t1n, x1] = ode45(@(t, q) odeFunctionRound1Negative(t, q, solutionRound1Positive), [t1p(end) 0], 0);
solutionRound1Negative = ode45(@(t, q) odeFunctionRound1Negative(t, q, solutionRound1Positive), [t1p(end) 0], 0);
x1_overflow_index = find(x1 < -2);
x1 = x1([1 x1_overflow_index(1)]);
t1n = t1n([1 x1_overflow_index(1)]);
figure
plot(t1n, x1, 'b-*')

function dqdt = odeFunctionRound1Positive(~, q)
    dqdt = pHq(-1, q);
end

function dxdt = odeFunctionRound1Negative(t, x, solution)
    q = deval(solution, linspace(t, t, 1));
    fprintf('t is %f, x is %f, we assume q is %f\n', t, x, q);
    dxdt = pHx(x, q);
end

function Hx = pHx(x, q)
    Hx = q - x - x.^3;
end

function Hq = pHq(x, q)
    Hq = q .* (3 * x.^2 + 1) - 8 * x.^3;
end
