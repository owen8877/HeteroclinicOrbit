param = getDefaultParam();

ySpan = linspace(0, 0.8, 101);
pySpan = linspace(-0.005, 0.005, 101);
hMat = zeros(size(pySpan, 2), size(ySpan, 2));

for i = 1:size(pySpan, 2)
    for j = 1:size(ySpan, 2)
        hMat(i, j) = Hrfunc(ySpan(j), pySpan(i), param);
    end
end

[y_, py_] = meshgrid(ySpan, pySpan);
% contour(y_, py_, hMat, 10, 'ShowText','on');
figure
surfc(y_, py_, hMat);

function Hr_ = Hrfunc(y, py, param)
    z = exp(py);
    b = param.b;
    Hr_ = (1/z - 1) * (y + (y + z / (b*(z-1)-1)) * (ft(y,param)-y*(1/z-1))/(gt(y,param)));
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

function Ft = ft(n, param)
    if nargin < 2
        error('Lacks param!');
    end
    K = param.K;
    Ft = f(n*K, param) / K;
end

function Gt = gt(n, param)
    if nargin < 2
        error('Lacks param!');
    end
    K = param.K;
    Gt = g(n*K, param) / K;
end