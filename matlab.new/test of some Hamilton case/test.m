clear;
clc;

xSolution = struct();
xSolution.t = [0 1];
xSolution.y = [-1 0];

C = 1;

xEnd = [];
errorT = [];
drawInterval = 1;

rounds = 20;

%% Constructing...

%% Iteration
for round = 1:rounds
    pSolution = myODESolver(@(s, p, x) pDfunc(s, p, x, C), [1 0], 1/2^round, xSolution, round);
    xSolution = myODESolver(@(s, x, p) xDfunc(s, x, p, C), [0 1], -1, pSolution, round);
    
    %%
    C = 0;
    for i = 1:(size(xSolution.y, 2) - 1)
        C = C + norm(xSolution.y(:, i) - xSolution.y(:, i+1));
    end
    fprintf('Round %d completed.\n', round);
    s_ = 0:0.01:1;
    x_ = interp1(xSolution.t', xSolution.y', s_, 'linear', 'extrap');
    p_ = interp1(pSolution.t', pSolution.y', s_, 'linear', 'extrap');
    errorT(round) = max(p_ - 2*(x_.^3 - x_));
    
    %%
    if round - fix(round / drawInterval) * drawInterval == 0
        figure
        subplot(1,2,1);
        plot(pSolution.t, pSolution.y);
        subplot(1,2,2);
        plot(xSolution.t, xSolution.y);
        title(['Round ' num2str(round)]);
%         plot(s_, p_ - 2*(x_.^3 - x_));
    end
end

figure
plot(1:1:rounds, log2(errorT));

function pD = pDfunc(s, p, x, C)
    f = p + x - x^3;
    g = p * (3*x^2 - 1);
    pD = g / norm(f) * C;
end

function xD = xDfunc(s, x, p, C)
    f = p + x - x^3;
    g = p * (3*x^2 - 1);
    xD = f / norm(f) * C;
end

function ySolution = myODESolver(yDfunc, tspan, y0, zSolution, round, varargin)
    %%%
    % The first row of `ySolution` is the array of time t
    %%%
    
    useZ = true;
    if nargin < 4
        zSolution = [];
        useZ = false;
    end
    
    tStart = tspan(1);
    tEnd = tspan(2);
    maxStep = abs(tEnd - tStart) / (2^7);
    minStep = abs(tEnd - tStart) / (2^23);
    adaptCoeff = 2^8;
    h = (tEnd - tStart) / (2^18);
    positiveDirection = tEnd > tStart;
    
    t = tStart;
    nextTimeQuit = false;
    
    index = 1;
    y = y0;
    
    if useZ
        tArray_zSolution = zSolution.t';
        zArray_zSolution = zSolution.y';
    end
    
    ySolution = struct();
    ySolution.t = [];
    ySolution.y = [];
    
    while true
        % adapt h
        if ~nextTimeQuit
            if abs(h) > maxStep
                h = sign(h) * maxStep;
            elseif abs(h) < minStep
                h = sign(h) * minStep;
            else
                if useZ
                    z = interp1(tArray_zSolution, zArray_zSolution, t, 'linear', 'extrap')';
                    yDAppr = yDfunc(t, y, z);
                else
                    yDAppr = yDfunc(t, y);
                end 

                h_ = norm(yDAppr) / adaptCoeff;
                while abs(h) < h_
                    h = h * 2;
                end
                while abs(h) > h_
                    h = h / 2;
                end
            end
        end
        
        ySolution.t(index) = t;
        ySolution.y(:, index) = y;
        
        if useZ
            z = interp1(tArray_zSolution, zArray_zSolution, t, 'linear', 'extrap')';
            z_ = interp1(tArray_zSolution, zArray_zSolution, t+h, 'linear', 'extrap')';
            y_ = y + h * yDfunc(t, y, z);
            y = y + h/2 * (yDfunc(t, y, z) + yDfunc(t+h, y_, z_));
        else
            y_ = y + h * yDfunc(t, y);
            y = y + h/2 * (yDfunc(t, y) + yDfunc(t+h, y_));
        end
        
        if nextTimeQuit
            break
        end
        
        if positiveDirection
            if t + h >= tEnd - h * 1e-3
                h = tEnd - t;
                nextTimeQuit = true;
            end
        else
            if t + h <= tEnd - h * 1e-3
                h = tEnd - t;
                nextTimeQuit = true;
            end
        end
        
        t = t + h;
        index = index + 1;
        clc;
        fprintf('round %d\t t %e\t h %e\n', round, t, h);
    end
end
