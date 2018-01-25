param = getDefaultParam();

timeLimit = 3;
resolution = 1e0;
dt = 1/resolution;
M = zeros(timeLimit*resolution, 1);
N = zeros(timeLimit*resolution, 1);
M(1) = 0;
N(1) = 0;

for i = 1:size(M, 1)-1
    m = M(i);
    n = N(i);
    F = f(n, param);
    G = g(n, param);
    M(i+1) = floor(max(m + (param.a*F/(F+G) - param.gamma*m)*1, 0));
    N(i+1) = floor(max(n + (param.gamma*param.b*m - n)*1, 0));
end

figure
subplot(2, 1, 1)
plot(1:timeLimit*resolution, M);
subplot(2, 1, 2)
plot(1:timeLimit*resolution, N);

function param = getDefaultParam()
%     param = struct();
%     a = 1600;
%     b = 2;
%     h = 2;
%     param.a = a;
%     param.b = b;
%     param.K = a*b;
%     param.k0min = a/50;
%     param.k0max = a;
%     param.k1min = a/100;
%     param.k1max = a/2;
%     param.gamma = 50;
%     param.h = h;
%     param.h1 = h;
%     param.h2 = h;
%     param.n50 = 2000;
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