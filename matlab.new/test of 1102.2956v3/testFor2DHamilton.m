clear
% clc
param = getDefaultParam();

L = @(y) y + (y - 1) * fTilde(y, param) / gTilde(y, param);
yZList = findAlmostAllZeroPoint(L, 0, 0.8, 0.005);

%% find the stable point in x-y-px-py space
yLeft = yZList(1);
yRight = yZList(2);
lStable = findZero(@(v) DH(v, param), [1/param.gamma; yLeft; 0; 0]);
lStable = lStable(1:2);
% lStable = [0; yLeft];
rStable = findZero(@(v) DH(v, param), [1/param.gamma; yRight; 0; 0]);
rStable = rStable(1:2);
% rStable = [0; yRight];
resolution = 1e5;
C = norm(lStable - rStable);
xySolution = zeros(2, resolution+1);
xySolution(1, :) = linspace(lStable(1), rStable(1), resolution+1);
xySolution(2, :) = linspace(lStable(2), rStable(2), resolution+1);
iteration = 0;
iterationLimit = 2;

while iteration < iterationLimit
    pSolution = naiveODESolver(@(s, p, xy) dpds(xy, p, param, C), ...
        [1 0], [0; 0], xySolution, resolution, @(p, xy) pcharDir(xy, p, param, C));
    xySolution = naiveODESolver(@(s, xy, p) dxyds(xy, p, param, C), ...
        [0 1], lStable, pSolution, resolution, @(xy, p) xycharDir(xy, p, param, C));
    
    C = 0;
    for i = 1:resolution
        C = C + norm(xySolution(:, i) - xySolution(:, i+1));
    end
    
    figure
    s1 = subplot(1, 2, 1); hold(s1, 'on');
    plot(xySolution(1, :), xySolution(2, :));
    plot(lStable(1), lStable(2), 'ro');
    plot(rStable(1), rStable(2), 'go');
    s2 = subplot(1, 2, 2); hold(s2, 'on');
    plot(pSolution(1, :), pSolution(2, :));
    plot(0, 0, 'ro');
    iteration = iteration + 1;
    clc
    fprintf('Iteration %d completed.\n', iteration);
end

%%
function charD = pcharDir(xy, p, param, C)
    H = @(v) Hfunc(v(1), v(2), v(3), v(4), param);
    hess = notOrdinaryHessian(H, [xy; p]);
    [charVec, D] = eig(hess);
    charD = - charVec(3:4, 4) / norm(charVec(1:2, 4)) * C;
end

function charD = xycharDir(xy, p, param, C)
    H = @(v) Hfunc(v(1), v(2), v(3), v(4), param);
    hess = notOrdinaryHessian(H, [xy; p]);
    [charVec, D] = eig(hess);
    charD = charVec(1:2, 3) / norm(charVec(1:2, 3)) * C;
end

function ySolution = naiveODESolver(yDfunc, tspan, y0, zSolution, resolution, charFunc)
    tStart = tspan(1);
    tEnd = tspan(2);
    h = (tEnd - tStart) / resolution;
    t = tStart;
    directionPos = tEnd > tStart;
    
    y = y0;
    ySolution = zeros(size(y0, 1), resolution+1);
    
    for index = 1:(resolution+1)
        if directionPos
            realIndex = index;
        else
            realIndex = resolution + 2 - index;
        end
        
        % record y
        ySolution(:, realIndex) = y;
        
        z = zSolution(:, realIndex);
        if index == 1
            yD = charFunc(y, z);
        else
            yk = y + h/2 * yDfunc(t, y, z);
            yD = yDfunc(t, yk, z);
        end
        y = y + h * yD;        
        t = t + h;
    end
end

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
        dh(i) = (Hfunc(vp(1), vp(2), vp(3), vp(4), param) - Hfunc(vm(1), vm(2), vm(3), vm(4), param)) / (2*d);
    end
end

function dxyds_ = dxyds(xy, p, param, C)
    dHdpx_ = dHdpx(xy(1), xy(2), p(1), p(2), param);
    dHdpy_ = dHdpy(xy(1), xy(2), p(1), p(2), param);
    dHdp = [dHdpx_; dHdpy_];
    H = Hfunc(xy(1), xy(2), p(1), p(2), param);
    dxyds_ = dHdp / sqrt(norm(dHdp)^2 + param.lambda * H) * C;
end

function dpds_ = dpds(xy, p, param, C)
    dHdpx_ = dHdpx(xy(1), xy(2), p(1), p(2), param);
    dHdpy_ = dHdpy(xy(1), xy(2), p(1), p(2), param);
    dHdp = [dHdpx_; dHdpy_];
    dHdx_ = dHdx(xy(1), xy(2), p(1), p(2), param);
    dHdy_ = dHdy(xy(1), xy(2), p(1), p(2), param);
    dHdxy = [dHdx_; dHdy_];
    H = Hfunc(xy(1), xy(2), p(1), p(2), param);
    dpds_ = - dHdxy / sqrt(norm(dHdp)^2 + param.lambda * H) * C;
end

function dHdx_ = dHdx(x, y, px, py, param)
    dx = 1e-7;
    dHdx_ = (Hfunc(x+dx, y, px, py, param) - Hfunc(x-dx, y, px, py, param)) / (2*dx);
%     gamma = param.gamma;
%     b = param.b;
%     pA = gamma*(exp(-px)-1)+gamma*b*(exp(py)-1);
%     A = y*(exp(-py)-1) + gamma*x*(exp(-px)-1) + gamma*b*x*(exp(py)-1);
%     hahaha = pA * (1+(fTilde(y, param)-2*A-(exp(px)-1)/b)/gTilde(y, param));
end

function dHdy_ = dHdy(x, y, px, py, param)
    dy = 1e-7;
    dHdy_ = (Hfunc(x, y+dy, px, py, param) - Hfunc(x, y-dy, px, py, param)) / (2*dy);
end

function dHdpx_ = dHdpx(x, y, px, py, param)
    dpx = 1e-7;
    dHdpx_ = (Hfunc(x, y, px+dpx, py, param) - Hfunc(x, y, px-dpx, py, param)) / (2*dpx);
end

function dHdpy_ = dHdpy(x, y, px, py, param)
    dpy = 1e-7;
    dHdpy_ = (Hfunc(x, y, px, py+dpy, param) - Hfunc(x, y, px, py-dpy, param)) / (2*dpy);
end

function H_ = Hfunc(x, y, px, py, param)
    gamma = param.gamma;
    b = param.b;
    A = y*(exp(-py)-1) + gamma*x*(exp(-px)-1) + gamma*b*x*(exp(py)-1);
    H_ = A + (A+(exp(px)-1)/b) * (fTilde(y, param)-A) / gTilde(y, param);
end

% function H__ = Hfunc(x, y, px, py, param)
%     z = exp(py);
%     b = param.b;
%     H__ = (1/z - 1) * (y + (y + z / (b*(z-1)-1)) * (fTilde(y,param)-y*(1/z-1))/(gTilde(y,param)));
% end

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