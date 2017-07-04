clear;
clc;

B = [[-1 4]; [-4 -1]];
% B = [[-1 0.5]; [-0.5 -1]];
p1 = [4; 1];
x0 = [0; 0];

xSolution = struct();
xSolution.t = [0 1];
xSolution.y = [[0; 1] [1; 0]];

C = sqrt(2);
% h = 1e-3;

xEnd = [];
drawInterval = 10;

%% Iteration
for round = 1:50
    pSolution = myODESolver(@(s, p, x) pDfunc(s, p, x, B, C), [0 1], p1, xSolution);
    xSolution = myODESolver(@(s, x, p) xDfunc(s, x, p, B, C), [1 0], [1/round; 0], pSolution);

    if round - fix(round / drawInterval) * drawInterval == 0
        figure
        title(['Round ' num2str(round)])
        subplot(1,2,1);
        plot(pSolution.y(1, :), pSolution.y(2, :));
        subplot(1,2,2);
        plot(xSolution.y(1, :), xSolution.y(2, :));
    end
    
    %%
    C = 0;
    for i = 1:(size(xSolution.y, 2) - 1)
        C = C + norm(xSolution.y(:, i) - xSolution.y(:, i+1));
    end
    fprintf('Round %d completed.\n', round);
    
    xEnd(:, round) = xSolution.y(:, end);
end

% for i = 1:399
%     n(i) = norm(xEnd(:, i) - xEnd(:, i+1));
% end
% 
% stem(log10(n));

function pD = pDfunc(s, p, x, B, C)
    b = - B * x;
    bM = norm(b);
    pD = C / bM * B' * p;
end

function xD = xDfunc(s, x, p, B, C)
    b = - B * x;
    bM = norm(b);
    xD = C / bM * (p + b);
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
    h = (tEnd - tStart) / 1e3;
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
            yD = yDfunc(t, y, z);
        else
            yD = yDfunc(t, y);
        end
        
        y = y + yD * h;
        
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
