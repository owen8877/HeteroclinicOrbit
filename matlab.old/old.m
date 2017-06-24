[X, Q] = meshgrid(linspace(-2, 1), linspace(-3, 1));

figure
streamslice(X, Q, pHx(X, Q), npHq(X, Q), 'arrow');
hold on;
plot(0, 0, 'r.', 'markersize', 20);
plot(-1, -2, 'b.', 'markersize', 20);
xlabel('x');
ylabel('q');
% hold off;

tMax = 5;
% interpolateSample = linspace(tMax, 0, 100);
tPositiveSpan = [0 tMax];
tNegativeSpan = [tMax 0];
epsilon = 1e-5;

opts = odeset('RelTol', 1e-2, 'AbsTol', 1e-4);

[t1p, x1] = ode45(@(t, x) pHx(x, 0), tPositiveSpan, -1, opts);
solutionForwards = ode45(@(t, x) pHx(x, 0), tPositiveSpan, -1);
% figure
% plot(t1p, x1, 'r-*')
% xlabel('t1p')
% ylabel('x1')
plot(x1, linspace(0, 0, length(x1)), 'r.')

[t1n, q1] = ode45(@(t, q) odeFuncn(t, q, solutionRound1Positive), tNegativeSpan, 0, opts);
solutionRound1Negative = ode45(@(t, q) odeFuncn(t, q, solutionRound1Positive), tNegativeSpan, 0);
% figure
% plot(t1n, q1, 'b-*')
% xlabel('t1n')
% ylabel('q1')
x1ForPlot = deval(solutionRound1Positive, t1n);
plot(x1ForPlot, q1, 'b.')

[t2p, x2] = ode45(@(t, x) odeFunc2p(t, x, solutionRound1Negative), tPositiveSpan, -1, opts);
solutionRound2Positive = ode45(@(t, x) odeFunc2p(t, x, solutionRound1Negative), tPositiveSpan, -1, opts);
% figure
% plot(t2p, x2, 'r-*')
% xlabel('t2p')
% ylabel('x2')
q1ForPlot = deval(solutionRound1Negative, t2p);
plot(x2, q1ForPlot, 'r.')

[t2n, q2] = ode45(@(t, q) odeFunc2n(t, q, solutionRound2Positive), tNegativeSpan, 0, opts);
solutionRound2Negative = ode45(@(t, q) odeFunc2n(t, q, solutionRound2Positive), tNegativeSpan, 0);
% figure
% plot(t1n, q1, 'b-*')
% xlabel('t1n')
% ylabel('q1')
x2ForPlot = deval(solutionRound2Positive, t2n);
plot(x2ForPlot, q2, 'b.')

function dxdt = odeFuncp(t, x, lastSolution)
    q = deval(lastSolution, t);
    % fprintf('t is %f, x is %f, we assume q is %f\n', t, x, q);
    dxdt = pHx(x, q);
end

function dqdt = odeFuncn(t, q, lastSolution)
    x = deval(lastSolution, t);
    % fprintf('t is %f, q is %f, we assume x is %f\n', t, q, x);
    dqdt = npHq(x, q);
end

function Hx = pHx(x, q)
    Hx = q - x - x.^3;
end

function Hq = npHq(x, q)
    Hq = q .* (3 * x.^2 + 1) - 8 * x.^3;
end