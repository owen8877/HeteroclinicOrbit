% [X, Q] = meshgrid(linspace(-2, 1), linspace(-3, 1));
% 
% figure
% streamslice(X, Q, pHx(X, Q), npHq(X, Q), 'arrow');
% hold on;
% plot(0, 0, 'r.', 'markersize', 20);
% plot(-1, -2, 'b.', 'markersize', 20);
% xlabel('x');
% ylabel('q');
% % hold off;

startPoint = [-1 -2];
endPoint = [0 0];

tMax = 20;
tPositiveSpan = [0 tMax];
tNegativeSpan = [tMax 0];

[t1p, x1] = ode45(@(t, x) pHx(x, 0), tPositiveSpan, startPoint(1));
solutionForwards = ode45(@(t, x) pHx(x, 0), tPositiveSpan, startPoint(1));
% figure
% plot(t1p, x1, 'r-*')
% xlabel('t1p')
% ylabel('x1')
% plot(x1, linspace(0, 0, length(x1)), 'r.')

lastx = x1;
lastq = [];

delta = [];

rounds = 1000;

for round = 1:rounds
    [tn, lastq] = ode45(@(t, q) odeFuncn(t, q, solutionForwards), tNegativeSpan, endPoint(2));
    solutionBackwards = ode45(@(t, q) odeFuncn(t, q, solutionForwards), tNegativeSpan, endPoint(2));
    % figure
    % plot(tn, lastq, 'b-*')
    % xlabel('tn')
    % ylabel('q')
    lastxForPlot = deval(solutionForwards, tn);
    delta(round) = max(lastq' - 2 * lastxForPlot .^ 3);
    % plot(lastxForPlot, lastq, 'b.')
    
    [tp, lastx] = ode45(@(t, x) odeFuncp(t, x, solutionBackwards), tPositiveSpan, startPoint(1));
    solutionForwards = ode45(@(t, x) odeFuncp(t, x, solutionBackwards), tPositiveSpan, startPoint(1));
    % figure
    % plot(tp, lastx, 'r-*')
    % xlabel('tp')
    % ylabel('x')
    lastqForPlot = deval(solutionBackwards, tp);
    % plot(lastx, lastqForPlot, 'r.')
end

semilogy(1:1:rounds, delta);

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
