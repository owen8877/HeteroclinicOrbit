function [lVal, rVal, Hqf, Hpf, Hf] = bioTestGen(default)
    if default
        getParameter(true);
        ly = 0.010526161790107;
        lx = 9.356588257872567e-06;
        lVal = [lx; ly; 0; 0];
        rx = 1.893653732450285e-04;
        ry = 0.213036044900658;
        rVal = [rx; ry; 0; 0];
    else
        getParameter(false);
        [lVal, rVal] = gridSearch();
    end
    Hqf = @Hqfunc;
    Hpf = @Hpfunc;
    Hf = @Hfunc;
end

function [lVal, rVal] = gridSearch()
    yWaitingList = [];
    phi = @(y) gTilde(y)*y - (1-y)*fTilde(y);
    
    figure
    hold on
    t = 0:1e-2:1;
    plot(t, gTilde(t).*t);
    plot(t, fTilde(t).*(1-t))
    legend('g*y', 'f*(1-y)')
    
    skipFlag = true;
    for y = linspace(0, 1, 21)
        p = phi(y);
        if skipFlag
            skipFlag = false;
            oldp = p;
            continue;
        end
        if sign(oldp) * sign(p) <= 0
            yWaitingList = [yWaitingList y];
        end
        oldp = p;
    end
    fprintf('Found %d candidates.\n', numel(yWaitingList))

    yCorrectList = [];
    for y = yWaitingList
        y1 = y;
        oldp = phi(y1);
        y2 = y + 1e-3;
        while true
            p = phi(y2);
            if norm(p) < 1e-16
                yCorrectList = [yCorrectList y2];
                break
            end
            
            y3 = y2 - p / (p-oldp) * (y2-y1);
            oldp = p;
            y1 = y2;
            y2 = y3;
        end
    end

    global bio_param
    gb = bio_param.gamma * bio_param.b;
    yC1 = yCorrectList(1);
    yC2 = yCorrectList(2);
    lVal = [yC1/gb; yC1; 0; 0]; rVal = [yC2/gb; yC2; 0; 0];
end

function dqHv = Hqfunc(q, p)
    h = 1e-8;
    dqHv = [ ...
        (Hfunc(q+[h; 0], p) - Hfunc(q-[h; 0], p)) / (2*h); ...
        (Hfunc(q+[0; h], p) - Hfunc(q-[0; h], p)) / (2*h) ...
        ];
end

function dpHv = Hpfunc(q, p)
    h = 1e-8;
    dpHv = [ ...
        (Hfunc(q, p+[h; 0]) - Hfunc(q, p-[h; 0])) / (2*h); ...
        (Hfunc(q, p+[0; h]) - Hfunc(q, p-[0; h])) / (2*h) ...
        ];
end

function H = Hfunc(q, p)
    global bio_param
    gamma = bio_param.gamma;
    b = bio_param.b;
    x = q(1); y = q(2);
    px = p(1); py = p(2);
    A = y*(exp(-py)-1) + gamma*x*(exp(-px)-1) + gamma*b*x*(exp(py)-1);
    H = A + (A+(exp(px)-1)/b) * (fTilde(y)-A) / gTilde(y);
end

function param = getParameter(default)
    param = struct();
    
    if default
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
    else
        a = 50;
        b = 40;
        h = 5;
        param.a = a;
        param.b = b;
        param.K = a*b;
        param.k0min = a/50;
        param.k0max = a;
        param.k1min = a/100;
        param.k1max = a;
        param.gamma = 5;
        param.h = h;
        param.h1 = h;
        param.h2 = h;
        param.n50 = 1000;
    end 
    
    global bio_param
    bio_param = param;
end

function F = f(n)
    global bio_param
    k0min = bio_param.k0min;
    k0max = bio_param.k0max;
    h1 = bio_param.h1;
    n50 = bio_param.n50;
    F = k0min + (k0max - k0min) * n.^h1 ./ (n50.^h1 + n.^h1);
end

function G = g(n)
    global bio_param
    k1min = bio_param.k1min;
    k1max = bio_param.k1max;
    h2 = bio_param.h2;
    n50 = bio_param.n50;
    G = k1max - (k1max - k1min) * n.^h2 ./ (n50.^h2 + n.^h2);
end

function Ft = fTilde(y)
    global bio_param
    K = bio_param.K;
    Ft = f(y*K) / K;
end

function Gt = gTilde(y)
    global bio_param
    K = bio_param.K;
    Gt = g(y*K) / K;
end