clear
% clc
param = getDefaultParam();

L = @(y) y + (y - 1) * fNorm(y, param) / gNorm(y, param);
yZList = findAlmostAllZeroPoint(L, 0, 0.8, 0.005);

%% we first find path from p1 to p2 above y-axis
yLeft = yZList(1);
yRight = yZList(2);
resolution = 1e3;
ySolution = linspace(yLeft, yRight, resolution+1);
C = yRight - yLeft;
iteration = 0;
iterationLimit = 100;
drawInterval = 0;

pErr = zeros(iterationLimit, 1);
while iteration < iterationLimit
    pSolution = naiveODESolver(@(s, p, y) dpds(y, p, param, C), [1 0], 0, ySolution, resolution, ...
        @(p, y) characterDirection(y, param));
    ySolution = naiveODESolver(@(s, y, p) dyds(y, p, param, C), [0 1], yLeft, pSolution, resolution, ...
        @(y, p) C);
    iteration = iteration + 1;
    
    C = 0;
    for i = 1:resolution
        C = C + norm(ySolution(:, i) - ySolution(:, i+1));
    end
    
    pPrecise = preciseSolution(ySolution, param);
    pErr(iteration) = max(abs(pSolution - pPrecise));
    
    if drawInterval > 0 && mod(iteration, drawInterval) == 0
        figure
        plot(ySolution, pSolution);
    end
    
    if true
        clc
        fprintf('Iteration %d / %d \n', iteration, iterationLimit);
    end
end

figure
plot(pErr)

figure
hold on
plot(ySolution, pSolution);
plot(ySolution, ySolution .* (ySolution .^ 2 - 1) * 2, 'r')

% [y, p] = meshgrid(0.2129:1e-5:0.2131, -1e-5:1e-6:1e-5);
% for i = 1:size(y, 1)
%     for j = 1:size(y, 2)
%         dy(i, j) = dpHr(y(i, j), p(i, j), param);
%         dp(i, j) = -dyHr(y(i, j), p(i, j), param);
%     end
% end
% quiver(y, p, dy, dp);

%%

function pPrecise = preciseSolution(ySpan, param)
    pPrecise = ySpan * 0;
    for i = 1:size(ySpan, 2)
        y = ySpan(i);
        b = param.b;
        A = (1+b*y)*(y+fNorm(y, param)) + b*y*gNorm(y, param);
        B = -y * (y*(1+2*b)+1+(1+b)*(fNorm(y, param)+gNorm(y, param)));
        C = (1+b) * y^2;
        pPrecise(i) = log((-B+sqrt(B^2-4*A*C))/(2*A));
    end
end

function xsList = findAlmostAllZeroPoint(f, a, b, scanStep)
    left = a;
    right = a + scanStep;
    
    xsList = [];
    while (right < b)
        if (f(left) * f(right) < 0)
            xsList = [xsList findZeroPoint(f, left, right)];
        end
        left = right;
        right = left + scanStep;
    end
end

function xs = findZeroPoint(f, a, b)
    dxlimit = 1e-5;
    if (a > b)
        error('The range boundary is reversed.')
    end
    
    while (b-a > dxlimit)
        xs = (a+b) / 2;
        if (f(a) * f(xs) > 0)
            a = xs;
        else
            b = xs;
        end
    end
end

function slope = characterDirection(y, param)
    slope = -2 * dpdyHrn(y, 0, param) / dpdpHrn(y, 0, param);
end

function dpdpHr_ = dpdpHrn(y, py, param)
    % use numeric method
    h = 1e-8;
    dpdpHr_ = (dpHr(y, py+h, param) - dpHr(y, py-h, param)) / (2*h);
end

function dpdyHr_ = dpdyHrn(y, py, param)
    % use numeric method
    h = 1e-8;
    dpdyHr_ = (dyHr(y, py+h, param) - dyHr(y, py-h, param)) / (2*h);
end

function dydyHr_ = dydyHrn(y, py, param)
    % use numeric method
    h = 1e-8;
    dydyHr_ = (dyHr(y+h, py, param) - dyHr(y-h, py, param)) / (2*h);
end

function dyds_ = dyds(y, py, param, C)
    dpHr_ = dpHr(y, py, param);
    H = Hrfunc(y, py, param);
    dyds_ = dpHr_ / sqrt(norm(dpHr_)^2 - param.lambda * H) * C;
end

function dpds_ = dpds(y, py, param, C)
    dpHr_ = dpHr(y, py, param);
    dyHr_ = dyHr(y, py, param);
    H = Hrfunc(y, py, param);
    dpds_ = - dyHr_ / sqrt(norm(dpHr_)^2 + param.lambda * H) * C;
end

function dpHr_ = dpHr(y, py, param)
    z = exp(py);
    b = param.b;
    dpHr_ = - y/z ...
        + (y + z/(b*(z-1)-1)) / gNorm(y, param) * (2*y/(z^2) - 2*y/z - fNorm(y, param)/z) ...
        - (1 - z) * (b+1) * (fNorm(y, param) - y * (1/z - 1)) / (b*(z-1)-1)^2 / gNorm(y, param);
end

function dyHr_ = dyHr(y, py, param)
    z = exp(py);
    b = param.b;
    K = param.K;
    h1 = param.h1;
    h2 = param.h2;
    n = K * y;
    n50 = param.n50;
    k0max = param.k0max;
    k0min = param.k0min;
    k1max = param.k1max;
    k1min = param.k1min;
    
    C1 = n50^h1;
    C2 = n50^h2;
    ftd = (k0max - k0min) * h1 * C1 * n^(h1-1) / (C1 + n^h1)^2;
    gtd = - (k1max - k1min) * h2 * C2 * n^(h2-1) / (C2 + n^h2)^2;
    ft = fNorm(y, param);
    gt = gNorm(y, param);
    dyHr_ = (1/z - 1) * ( ...
        1 + (ft - y * (1/z-1)) / gt ...
        + (y+z/(b*(z-1)-1)) * ((ftd-(1/z-1))*gt - (ft-y*(1/z-1))*gtd) / gt^2 );
end

function Hr_ = Hrfunc(y, py, param)
    z = exp(py);
    b = param.b;
    Hr_ = (1/z - 1) * (y + (y + z / (b*(z-1)-1)) * (fNorm(y,param)-y*(1/z-1))/(gNorm(y,param)));
end

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
    param.lambda = 1e-2;
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

function Ft = fNorm(y, param)
    if nargin < 2
        error('Lacks param!');
    end
    K = param.K;
    Ft = f(y*K, param) / K;
end

function Gt = gNorm(y, param)
    if nargin < 2
        error('Lacks param!');
    end
    K = param.K;
    Gt = g(y*K, param) / K;
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
            yD = yDfunc(t, y, z);
        end
        y = y + h * yD;        
        t = t + h;
    end
end