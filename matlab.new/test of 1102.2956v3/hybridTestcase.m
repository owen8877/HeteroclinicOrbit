clear
clc
param = getDefaultParam();

L = @(y) y + (y - 1) * fTilde(y, param) / gTilde(y, param);
yZList = findAlmostAllZeroPoint(L, 0, 0.8, 0.005);

%% find the stable point in x-y-px-py space
yLeft = yZList(1);
yRight = yZList(2);
lStable = findZero(@(v) DH(v, param), [1/param.gamma; yLeft; 0; 0]);
lStable = lStable(1:2);
rStable = findZero(@(v) DH(v, param), [1/param.gamma; yRight; 0; 0]);
rStable = rStable(1:2);
resolution = 1e4;
C = norm(lStable - rStable);
iteration = 0;
iterationLimit = 1;

return

lHess = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4), param), [lStable; 0; 0]);
rHess = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4), param), [rStable; 0; 0]);

%%
[V, D] = eig(lHess);

%%
Tspan = 10;
scale = 1 / resolution;
% delta = rand(4,1);
% delta(4) = -abs(delta(4));
% delta(3) = -abs(delta(3));
% delta = delta / norm(delta);
delta = V(:, 4);

% solution = simpleSymplecticODESolver(...
%     @(~, v, ~) dH(v(1:2), v(3:4), param), [0, Tspan], ...
%     [lStable; 0; 0] + delta * scale, [], resolution);
solution = simpleSymplecticSearch(...
    @(~, v, ~) dH(v(1:2), v(3:4), param), ...
    [lStable; 0; 0] + delta * scale, [rStable; 0; 0], ...
    1 / resolution, 0);

figure
plotHelper(subplot(2, 2, 1), solution([1, 3], :), 'x', 'px', [[lStable(1) rStable(1)]; [0 0]]);
plotHelper(subplot(2, 2, 2), solution([2, 4], :), 'y', 'py', [[lStable(2) rStable(2)]; [0 0]]);
plotHelper(subplot(2, 2, 3), solution([1, 2], :), 'x', 'y', [lStable rStable]);
plotHelper(subplot(2, 2, 4), solution([3, 4], :), 'px', 'py', [0; 0]);

H = zeros(1, resolution+1);
for i = 1:resolution+1
    H(i) = Hfunc(solution(1:2, i), solution(3:4, i), param);
end

figure
plot(abs(H))

%%
function pointList = findAlmostAllZeroPoint(f, lRange, rRange, step)
    pointList = [];
    leftSign = f(lRange) > 0;
    while lRange < rRange
        lRange = lRange + step;
        rightSign = f(lRange) > 0;
        if xor(leftSign, rightSign)
            pointList = [pointList findZero(f, lRange)];
        end
        leftSign = rightSign;
    end
end

function dh = DH(v, param)
    dh = v * 0;
    d = 1e-7;
    for i = 1:size(dh, 1)
        vp = v; vp(i) = v(i) + d;
        vm = v; vm(i) = v(i) - d;
        dh(i) = (Hfunc(vp(1:2), vp(3:4), param) - Hfunc(vm(1:2), vm(3:4), param)) / (2*d);
    end
end

% function dxdsv = dxds(x, p, param, C)
%     dpHv = dpHn(x, p);
%     H = Hfunc(x, p);
%     dxdsv = dpHv / sqrt(norm(dpHv, 2)^2 - param.lambda * H) * C;
% end
% 
% function dpdsv = dpds(x, p, param, C)
%     dpHv = dpHn(x, p);
%     dxHv = dxHn(x, p);
%     H = Hfunc(x, p);
%     dpdsv = - dxHv / sqrt(norm(dpHv, 2)^2 - param.lambda * H) * C;
% end

function dHv = dH(x, p, param)
    dxH = dxHnum(x, p, param);
    dpH = dpHnum(x, p, param);
    dHv = [dpH; -dxH];
    dHv = dHv / norm(dHv);
end

function dHv = dHsimple(x, p, param)
    dxH = dxHnum(x, p, param);
    dpH = dpHnum(x, p, param);
    dHv = [dpH; -dxH];
end

function dxxHv = dxxHnum(x, p, param)
    h = 1e-10;
    dxxHv = [ ...
        (dxHnum(x+[h; 0], p, param) - dxHnum(x-[h; 0], p, param)) / (2*h) ...
        (dxHnum(x+[0; h], p, param) - dxHnum(x-[0; h], p, param)) / (2*h) ...
        ];
end

function dpxHv = dpxHnum(x, p, param)
    h = 1e-10;
    dpxHv = [ ...
        (dpHnum(x+[h; 0], p, param) - dpHnum(x-[h; 0], p, param)) / (2*h) ...
        (dpHnum(x+[0; h], p, param) - dpHnum(x-[0; h], p, param)) / (2*h) ...
        ];
end

function dxpHv = dxpHnum(x, p, param)
    dxpHv = dpxHnum(x, p, param)';
end

function dppHv = dppHnum(x, p, param)
    h = 1e-10;
    dppHv = [ ...
        (dpHnum(x, p+[h; 0], param) - dpHnum(x, p-[h; 0], param)) / (2*h) ...
        (dpHnum(x, p+[0; h], param) - dpHnum(x, p-[0; h], param)) / (2*h) ...
        ];
end

function dxHv = dxHnum(x, p, param)
    h = 1e-10;
    dxHv = [ ...
        (Hfunc(x+[h; 0], p, param) - Hfunc(x-[h; 0], p, param)) / (2*h); ...
        (Hfunc(x+[0; h], p, param) - Hfunc(x-[0; h], p, param)) / (2*h) ...
        ];
end

function dpHv = dpHnum(x, p, param)
    h = 1e-10;
    dpHv = [ ...
        (Hfunc(x, p+[h; 0], param) - Hfunc(x, p-[h; 0], param)) / (2*h); ...
        (Hfunc(x, p+[0; h], param) - Hfunc(x, p-[0; h], param)) / (2*h) ...
        ];
end

function H = Hfunc(xy, p, param)
    gamma = param.gamma;
    b = param.b;
    x = xy(1); y = xy(2);
    px = p(1); py = p(2);
    A = y*(exp(-py)-1) + gamma*x*(exp(-px)-1) + gamma*b*x*(exp(py)-1);
    H = A + (A+(exp(px)-1)/b) * (fTilde(y, param)-A) / gTilde(y, param);
end

%%
function param = getDefaultParam()
    param = struct();
    a = 320/3;
    b = 22.5;
    h = 2;
    param.a = a;
    param.b = b;
    param.K = a*b;
    param.k0min = a/100;
    param.k0max = a;
    param.k1min = a/100;
    param.k1max = a;
    param.gamma = 50;
    param.h = h;
    param.h1 = h;
    param.h2 = h;
    param.n50 = 1000;
    param.lambda = 2e-3;
end

%%
function F = f(n, param)
    if nargin < 2
        error('Lacks param!');
    end
    k0max = param.k0max;
    k0min = param.k0min;
    h1 = param.h1;
    n50 = param.n50;
    F = k0min + (k0max - k0min) * n^h1 / (n50^h1 + n^h1);
end

function G = g(n, param)
    if nargin < 2
        error('Lacks param!');
    end
    k1max = param.k1max;
    k1min = param.k1min;
    h2 = param.h2;
    n50 = param.n50;
    G = k1max - (k1max - k1min) * n^h2 / (n50^h2 + n^h2);
end

function Ft = fTilde(y, param)
    if nargin < 2
        error('Lacks param!');
    end
    K = param.K;
    Ft = f(y*K, param) / K;
end

function Gt = gTilde(y, param)
    if nargin < 2
        error('Lacks param!');
    end
    K = param.K;
    Gt = g(y*K, param) / K;
end