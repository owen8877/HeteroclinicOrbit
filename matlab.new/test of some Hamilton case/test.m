clear;
clc;

% p1 = [4; 1];
% x0 = [0; 0];

xSolution = struct();
xSolution.t = [0 1];
xSolution.y = [0 -1];

C = 1;
% h = 1e-3;

xEnd = [];
errorT = [];
drawInterval = 1;

%% Iteration
for round = 1:10
    pSolution = myODESolver(@(s, p, x) pDfunc(s, p, x, C), [0 1], -1/2^round, xSolution);
    xSolution = myODESolver(@(s, x, p) xDfunc(s, x, p, C), [1 0], -1, pSolution);

    if round - fix(round / drawInterval) * drawInterval == 0
        figure
        subplot(1,2,1);
        plot(pSolution.t, pSolution.y);
        subplot(1,2,2);
        plot(xSolution.t, xSolution.y);
        title(['Round ' num2str(round)]);
    end
    
    %%
    C = 0;
    for i = 1:(size(xSolution.y, 2) - 1)
        C = C + norm(xSolution.y(:, i) - xSolution.y(:, i+1));
    end
    fprintf('Round %d completed.\n', round);
    s_ = 0:0.01:1;
    x_ = interp1(xSolution.t', xSolution.y', s_, 'linear', 'extrap');
    p_ = interp1(pSolution.t', pSolution.y', s_, 'linear', 'extrap');
    errorT(round) = max(p_);
end

figure
plot(1:1:10, log10(errorT));

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

function ySolution = myODESolver(yDfunc, tspan, y0, zSolution, varargin)
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
    h = (tEnd - tStart) / (2^10);
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
    end
end
