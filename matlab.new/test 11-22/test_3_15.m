clc; clear

% init
q1Stable = -1/sqrt(2);
q1Saddle = 0;
% % for bio cases
% ly = 0.010526161790107;
% lx = 9.356588257872567e-06;
% lVal = [lx; ly; 0; 0];
% rx = 1.893653732450285e-04;
% ry = 0.213036044900658;
% rVal = [rx; ry; 0; 0];
% % for prob cases
% lVal = [q1Stable; 0; 0; 0];
% rVal = [q1Saddle; 0; 0; 0];
% for self-made cases
lVal = [-1; -1; 0; 0];
rVal = [0; 0; 0; 0];

resolution = 3000;

lHessian = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4)), lVal);
% lHessian = nddH(lVal);
[lV, lD] = eigs(lHessian, size(lHessian, 1), 'lr');
rHessian = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4)), rVal);
% rHessian = nddH(rVal);
[rV, rD] = eigs(rHessian, size(rHessian, 1), 'lr');
leftSide = false;
gridSearch = false;

%% here the right anchor has two stable manifold, which can take as tests.
% lAnchor = lVal + real(lV(:, 4)) / resolution;
% solution = simpleSymplecticSearch(@(~, v, ~) dHn(v(1:2), v(3:4)), ...
%     lAnchor, rVal, 1/resolution, 0, rV(:, end));
% % richPlotHelper(0, solution, lVal, rVal, rVal, @Hfunc, @dH);
% 
% fai = atan(2.8);
% rAnchor = rVal + (cos(fai) * rV(:, 1) + sin(fai) * rV(:, 2)) / resolution;
% initialSolution = pulling(solution, 0.05, 3, rAnchor);
% % richPlotHelper(0, initialSolution, lVal, rVal, rAnchor, @Hfunc, @dH);
% 
% T = (size(initialSolution, 2) - 1) / resolution;
% q1Solution = initialSolution;
% q1Solution(1:2, :) = simpleSymplecticODESolver(@(t, q, p) dHdqn(q, p), ...
%     [T 0], rAnchor(1:2), initialSolution(3:4, :), size(initialSolution, 2) - 1);
% richPlotHelper(0, q1Solution, lVal, rVal, rAnchor, @Hfunc, @dH);

%% The following code gives a set of possible initial angle between eigen directions
fais = linspace(0, 2*pi, 100)';
fais = fais(2:end);
distance = [];
notice = [];
if gridSearch
    for fai = fais'
        [thisDistance, ~] = simplexSearchWrapper(fai, leftSide, lVal, rVal, ...
            lV, rV, lD, rD, @dHn, resolution);
        distance = [distance; thisDistance];
        fprintf('Searching %.6e with distance %.4e\n', fai, thisDistance);
        if thisDistance < 1
            notice = [notice; fai];
            fprintf('Found one %.6e\n', fai);
        end
    end
    figure
    plot(fais, distance);
end

%% The following code searches for the minimal value
% initial = notice(1);
% initial = 0.4649557127; % for self-made cases on the right, L=3
initial = 0.3769; % for self-made cases on the right, L=5
% initial = 1.55823; % for prob cases
% initial = 4; % for bio cases
% initial = 2; % for self-made cases on the left, L=3
step = pi/1000;

[initialDis, ~] = simplexSearchWrapper(initial, leftSide, lVal, rVal, lV, rV, ...
    lD, rD, @dHn, resolution);
fprintf('Init dis %.4e\n', initialDis);

for leftMulti = 1:20
    fai = initial - (leftMulti+1) * step;
    [dis, ~] = simplexSearchWrapper(fai, leftSide, lVal, rVal, lV, rV, ...
        lD, rD, @dHn, resolution);
    fprintf('LeftMulti %d distance %.4e\n', leftMulti, dis);
    if dis > initialDis
        break
    end
end
for i = 1:10
    left = initial - (leftMulti+2^(-i)) * step;
    [dis, ~] = simplexSearchWrapper(left, leftSide, lVal, rVal, lV, rV, ...
        lD, rD, @dHn, resolution);
    fprintf('Lefti %d distance %.4e\n', i, dis);
    if dis < initialDis
        break
    end
end

for rightMulti = 1:20
    fai = initial + (rightMulti+1) * step;
    [dis, ~] = simplexSearchWrapper(fai, leftSide, lVal, rVal, lV, rV, ...
        lD, rD, @dHn, resolution);
    fprintf('RightMulti %d distance %.4e\n', rightMulti, dis);
    if dis > initialDis
        break
    end
end
for i = 1:10
    right = initial + (rightMulti+2^(-i)) * step;
    [dis, ~] = simplexSearchWrapper(right, leftSide, lVal, rVal, lV, rV, ...
        lD, rD, @dHn, resolution);
    fprintf('Righti %d distance %.4e\n', i, dis);
    if dis < initialDis
        break
    end
end

leftValueKnown = false; rightValueKnown = false;
for step = 1:20
    midleft = right - (right-left) * 0.618;
    midright = left + (right-left) * 0.618;
    if ~leftValueKnown
        [midleftValue, solution] = simplexSearchWrapper(midleft, leftSide, ...
            lVal, rVal, lV, rV, lD, rD, @dHn, resolution);
    end
    
    if ~rightValueKnown
        [midrightValue, solution] = simplexSearchWrapper(midright, leftSide, ...
            lVal, rVal, lV, rV, lD, rD, @dHn, resolution);
    end
    
    fprintf('%.6e, %.6e - %.4e, %.4e\n', left, right, midleftValue, midrightValue);
        
    if midleftValue < midrightValue
        right = midright;
        midrightValue = midleftValue;
        rightValueKnown = true;
        leftValueKnown = false;
    else
        left = midleft;
        midleftValue = midrightValue;
        rightValueKnown = false;
        leftValueKnown = true;
    end
end
richPlotHelper(0, solution, lVal, rVal, lVal, @Hfunc, @dH);

%%
function [distance, solution] = simplexSearchWrapper(fai, leftSide, lVal, rVal, lV, rV, lD, rD, dHn, resolution)
    if ~leftSide
        rAnchor = rVal + (rV(:, 1)*sin(fai) + rV(:, 2)*cos(fai)) / resolution;
        solution = simpleSymplecticSearch(@(~, v, ~) dHn(v(1:2), v(3:4)), ...
            rAnchor, lVal, -1/resolution, 0, real(lV(:, end)), struct('output', false));
        solution = fliplr(solution);
        % richPlotHelper(0, solution, lVal, rVal, lVal, @Hfunc, @dH);
        distance = norm(lVal - solution(:, 1));
    else
        if isreal(lD(3, 3))
            delta = lV(:, 4)*sin(fai) + lV(:, 3)*cos(fai);
        else
            delta = lV(:, 4)*sin(fai) + imag(lV(:, 3))*cos(fai);
            delta = delta / norm(delta);
        end
        lAnchor = lVal + delta / resolution;
        solution = simpleSymplecticSearch(@(~, v, ~) dHn(v(1:2), v(3:4)), ...
            lAnchor, rVal, 1/resolution, 2, real(rV(:, 1)), struct('output', false));
        distance = norm(rVal - solution(:, end));
    end
end

%%
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
% function nddHp = nddH(v)
%     q1 = v(1); q2 = v(2); p1 = v(3); p2 = v(4);
%     error('blabla');
%     nddHp = [0 1 0 0; ddf(q1)*q2 df(q1) 0 0.5; ...
%         (-p2*q2*dddf(q1)+2*ddf(q1)*f(q1)+2*(df(q1))^2-2*ddf(q1)*q2) ...
%         (-p2*ddf(q1)-2*df(q1)) 0 (-ddf(q1)*q2); ...
%         (-p2*ddf(q1)-2*df(q1)) 2 -1 -df(q1)];
% end

function Hpv = Hpfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    Hpv = [p1+q1-q1^3; p2+q2-q2^3];
end

function Hqv = Hqfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    L = 5;
    Hqv = [p1*(1-3*q1^2)+2*((q2^2-1)*q1^2-L*(q1^2-1)*q2^2)*(2*q1*(q2^2-1)-L*2*q1*q2^2); ...
        p2*(1-3*q2^2)+2*((q2^2-1)*q1^2-L*(q1^2-1)*q2^2)*(2*q2*q1^2-L*2*q2*(q1^2-1))];
end

function H = Hfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    L = 5;
    H = (p1^2+p2^2)/2 + p1*(q1-q1^3) + p2*(q2-q2^3) + ...
        ((q2^2-1)*q1^2-L*(q1^2-1)*q2^2)^2;
end

%% prob H func
% needs checking
% function nddHp = nddH(v)
%     q1 = v(1); q2 = v(2); p1 = v(3); p2 = v(4);
%     nddHp = [0 1 0 0; ddf(q1)*q2 df(q1) 0 0.5; ...
%         (-p2*q2*dddf(q1)+2*ddf(q1)*f(q1)+2*(df(q1))^2-2*ddf(q1)*q2) ...
%         (-p2*ddf(q1)-2*df(q1)) 0 (-ddf(q1)*q2); ...
%         (-p2*ddf(q1)-2*df(q1)) 2 -1 -df(q1)];
% end
% 
% function Hpv = Hpfunc(q, p)
%     p1 = p(1); p2 = p(2);
%     q1 = q(1); q2 = q(2);
%     Hpv = [q2; p2/2+df(q1)*q2];
% end
% 
% function Hqv = Hqfunc(q, p)
%     p1 = p(1); p2 = p(2);
%     q1 = q(1); q2 = q(2);
%     Hqv = [p2*q2*ddf(q1)+2*(q2-f(q1))*df(q1); p1+p2*df(q1)-2*(q2-f(q1))];
% end
% 
% function H = Hfunc(q, p)
%     p1 = p(1); p2 = p(2);
%     q1 = q(1); q2 = q(2);
%     H = p1*q2 + p2*(p2/4+df(q1)*q2) - (q2-f(q1))^2;
% end
% 
% function fv = f(x)
%     % fv = -(5*x^4 + 8*x^3 + 3*x^2);
%     % fv = 12*x^3 + 12*x^2;
%     % fv = 2*x^3 - x/2;
%     fv = 4*x^3 - 2*x;
% end
% 
% function dfv = df(x)
%     % dfv = -(20*x^3 + 24*x^2 + 6*x);
%     % dfv = 36*x^2 + 24;
%     % dfv = 6*x^2 - 1/2;
%     dfv = 12*x^2 - 2;
% end
% 
% function ddfv = ddf(x)
%     % ddfv = -(60*x^2 + 48*x + 6);
%     % ddfv = 72*x + 24;
%     % ddfv = 12*x;
%     ddfv = 24*x;
% end
% 
% function dddfv = dddf(x)
%     % dddfv = -(120*x + 48);
%     % dddfv = 72;
%     % dddfv = 12;
%     dddfv = 24;
% end

%% bio H func
% function dqHv = Hqfunc(q, p)
%     h = 1e-8;
%     dqHv = [ ...
%         (Hfunc(q+[h; 0], p) - Hfunc(q-[h; 0], p)) / (2*h); ...
%         (Hfunc(q+[0; h], p) - Hfunc(q-[0; h], p)) / (2*h) ...
%         ];
% end
% 
% function dpHv = Hpfunc(q, p)
%     h = 1e-8;
%     dpHv = [ ...
%         (Hfunc(q, p+[h; 0]) - Hfunc(q, p-[h; 0])) / (2*h); ...
%         (Hfunc(q, p+[0; h]) - Hfunc(q, p-[0; h])) / (2*h) ...
%         ];
% end
% 
% function H = Hfunc(q, p)
%     gamma = 50;
%     b = 22.5;
%     x = q(1); y = q(2);
%     px = p(1); py = p(2);
%     A = y*(exp(-py)-1) + gamma*x*(exp(-px)-1) + gamma*b*x*(exp(py)-1);
%     H = A + (A+(exp(px)-1)/b) * (fTilde(y)-A) / gTilde(y);
% end
% 
% function param = getDefaultParam()
%     param = struct();
%     a = 320/3;
%     b = 22.5;
%     h = 2;
%     param.a = a;
%     param.b = b;
%     param.K = a*b;
%     param.k0min = a/100;
%     param.k0max = a;
%     param.k1min = a/100;
%     param.k1max = a;
%     param.gamma = 50;
%     param.h = h;
%     param.h1 = h;
%     param.h2 = h;
%     param.n50 = 1000;
%     param.lambda = 1e-2;
% end
% 
% function F = f(n)
%     a = 320/3;
%     k0min = a/100;
%     k0max = a;
%     h = 2;
%     h1 = h;
%     n50 = 1000;
%     F = k0min + (k0max - k0min) * n^h1 / (n50^h1 + n^h1);
% end
% 
% function G = g(n)
%     a = 320/3;
%     k1min = a/100;
%     k1max = a;
%     h = 2;
%     h2 = h;
%     n50 = 1000;
%     G = k1max - (k1max - k1min) * n^h2 / (n50^h2 + n^h2);
% end
% 
% function Ft = fTilde(y)
%     a = 320/3;
%     b = 22.5;
%     K = a*b;
%     Ft = f(y*K) / K;
% end
% 
% function Gt = gTilde(y)
%     a = 320/3;
%     b = 22.5;
%     K = a*b;
%     Gt = g(y*K) / K;
% end