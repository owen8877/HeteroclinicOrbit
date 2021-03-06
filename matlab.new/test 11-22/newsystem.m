clc; clear;

% init
q1Stable = -1/sqrt(2);
q1Saddle = 0;
lVal = [q1Stable; 0; 0; 0];
rVal = [q1Saddle; 0; 0; 0];
rVal2 = rVal;
resolution = 1e3;

% lHessian = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4)), lVal);
lHessian = nddH(lVal);
[lV, lD] = eig(lHessian);
% rHessian = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4)), rVal);
rHessian = nddH(rVal);
[rV, rD] = eig(rHessian);

rAnchor = rVal + real(rV(:, 3)) / resolution;
lAnchor = lVal - real(lV(:, 3)) / resolution;

%% Search part
V = lV; V(:, [1 3]) = lV(:, [3 1]);
solution = simpleSymplecticSearch(@(~, v, ~) dHn(v(1:2), v(3:4)), ...
    rAnchor, lVal, -1/resolution, 1, V);

% newSolution = solution;
initialSolution = pulling(solution, 0.5, 3, lAnchor);
initialSolution = fliplr(initialSolution);

richPlotHelper(0, initialSolution, lVal, rVal, rVal2, @Hfunc, @dH);

%%
function dHv = dH(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
end

function dHnv = dHn(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
    dHnv = dHv / norm(dHv);
end

function dHdpnv = dHdpn(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
    dHnv = dHv / norm(dHv);
    dHdpnv = dHnv(3:4);
end

function dHdqnv = dHdqn(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
    dHnv = dHv / norm(dHv);
    dHdqnv = dHnv(1:2);
end

function nddHp = nddH(v)
    q1 = v(1); q2 = v(2); p1 = v(3); p2 = v(4);
    nddHp = [0 1 0 0; ddf(q1)*q2 df(q1) 0 0.5; ...
        (-p2*q2*dddf(q1)+2*ddf(q1)*f(q1)+2*(df(q1))^2-2*ddf(q1)*q2) (-p2*ddf(q1)-2*df(q1)) 0 (-ddf(q1)*q2); ...
        (-p2*ddf(q1)-2*df(q1)) 2 -1 -df(q1)];
end

function Hpv = Hpfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    Hpv = [q2; p2/2+df(q1)*q2];
end

function Hqv = Hqfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    Hqv = [p2*q2*ddf(q1)+2*(q2-f(q1))*df(q1); p1+p2*df(q1)-2*(q2-f(q1))];
end

function H = Hfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    H = p1*q2 + p2*(p2/4+df(q1)*q2) - (q2-f(q1))^2;
end

function fv = f(x)
    % fv = -(5*x^4 + 8*x^3 + 3*x^2);
    % fv = 12*x^3 + 12*x^2;
    % fv = 2*x^3 - x/2;
    fv = 4*x^3 - 2*x;
end

function dfv = df(x)
    % dfv = -(20*x^3 + 24*x^2 + 6*x);
    % dfv = 36*x^2 + 24;
    % dfv = 6*x^2 - 1/2;
    dfv = 12*x^2 - 2;
end

function ddfv = ddf(x)
    % ddfv = -(60*x^2 + 48*x + 6);
    % ddfv = 72*x + 24;
    % ddfv = 12*x;
    ddfv = 24*x;
end

function dddfv = dddf(x)
    % dddfv = -(120*x + 48);
    % dddfv = 72;
    % dddfv = 12;
    dddfv = 24;
end