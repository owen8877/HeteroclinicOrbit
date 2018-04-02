clc; clear

load initialSolution.mat

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

%%
solution = initialSolution;
richPlotHelper(0, solution, lVal, rVal, rVal2, @Hfunc, @dH);

for itr = 1:1000
    for i = 2:size(solution, 2)-1
        qp = solution(:, i);
        qpprevious = solution(:, i-1);
        qpnext = solution(:, i+1);
        updateD = - 0.1 * gradient(qp, qpprevious, qpnext);
        solution(:, i) = qp + updateD;
    end
    if mod(itr, 10) == 0
%     	richPlotHelper(itr, solution, lVal, rVal, rVal2, @Hfunc, @dH);
    end
    if mod(itr, 10) == 0
        % clc
        fprintf('Itr\t%d\n', itr);
    	richReportHelper(itr, solution, lVal, rVal, rVal2, @Hfunc, @dH);
    end
end
richPlotHelper(itr, solution, lVal, rVal, rVal2, @Hfunc, @dH);

%%
function g = gradient(v, vprevious, vnext)
    lambda = 5;
    gH = nablaH(v(1:2), v(3:4));
    
    hz = dH(vprevious(1:2), vprevious(3:4));
    deltaz = v - vprevious;
    w = hz*norm(deltaz) - deltaz*norm(hz);
    nablahz = nddH(vprevious);
%     nablaw = nablahz*norm(deltaz) + hz*deltaz'/norm(deltaz) ...
%         - eye(4)*norm(hz) - deltaz * (hz' * nablahz) / norm(hz);
    nablaw = hz*deltaz'/norm(deltaz) - eye(4)*norm(hz);

    hz_ = dH(vnext(1:2), vnext(3:4));
    deltaz_ = vnext - v;
    w_ = hz_*norm(deltaz_) - deltaz_*norm(hz_);
    nablahz_ = nddH(v);
%     nablaw_ = nablahz_*norm(deltaz_) - hz_*deltaz_'/norm(deltaz_) ...
%         + eye(4)*norm(hz_) - deltaz_ * (hz_' * nablahz_) / norm(hz_);
    nablaw_ = 0 - hz_*deltaz_'/norm(deltaz_) + eye(4)*norm(hz_);

    % g = 2*gH*Hfunc(v(1:2), v(3:4)) + lambda * 2 * (w'*nablaw + w_'*nablaw_)';
    g = lambda * 2 * (w'*nablaw + w_'*nablaw_)';
    % g = 2*gH*Hfunc(v(1:2), v(3:4)) + lambda * 2 * (w'*nablaw)';

%     dHz = dH(vprevious(1:2), vprevious(3:4));
%     dz = v - vprevious;
%     Lz = dHz - normalize(dz)*norm(dHz);
%     r = norm(dz);
%     dLz = 0 - (eye(4)/r-dz*dz'/(r^3)) * norm(dHz);
%     g = 2*gH*Hfunc(v(1:2), v(3:4)) + lambda*2*(dLz*Lz);
end

function nv = normalize(v)
    nv = v / norm(v);
end

function dHv = dH(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
end

function dHv = nablaH(q, p)
    dHv = [Hqfunc(q, p); Hpfunc(q, p)];
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

% needs checking
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