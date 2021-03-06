clear;
clc;

B = [[-1 4]; [-4 -1]];
p0 = [4; 1];
x1 = [0; 0];

drawInterval = 1;
iterationLimit = 10;

resolution = 10^3;
criterion = 1e-7;
iteration = 0;
p1Collection = [];
previousP1 = [];

sArray = 0:1/resolution:1;
% init xSolution
xSolution = [1-sArray; sArray];
C = sqrt(2);

%% Iteration
while true
    pSolution = naiveODESolver(@(s, p, x) pDfunc(s, p, x, B, C), [0 1], p0, xSolution, resolution);
    xSolution = naiveODESolver(@(s, x, p) xDfunc(s, x, p, B, C), [1 0], x1, pSolution, resolution);

    iteration = iteration + 1;
    if iteration - fix(iteration / drawInterval) * drawInterval == 0
        figure;
        subplot(1,2,1);
        plot(pSolution(1, :), pSolution(2, :));
        title(['pSolution, iteration ' num2str(iteration)])

        subplot(1,2,2);
        plot(xSolution(1, :), xSolution(2, :));
        title(['xSolution, iteration ' num2str(iteration)])
    end

    %% update C
    C = 0;
    for i = 1:(size(xSolution, 2) - 1)
        C = C + norm(xSolution(:, i) - xSolution(:, i+1));
    end

    %% output iteration info
    clc
    fprintf('(Resolution %.0e) Iteration %d completed.\n', resolution, iteration);

    if iteration > iterationLimit
        break
    end

    %% watch if p(1) converges
    if iteration ~= 1
        p1 = pSolution(:, end);
        endDistance = norm(p1 - previousP1);
        if endDistance < criterion
            break
        end
        fprintf('probable error in p1: %e', endDistance);
    end
    previousP1 = pSolution(:, end);
    p1Collection(:, iteration) = previousP1;
end

function pD = pDfunc(~, p, x, B, C)
    b = - B * x;
    bM = norm(b);
    pD = C / bM * B' * p;
end

function xD = xDfunc(~, x, p, B, C)
    b = - B * x;
    bM = norm(b);
    xD = C / bM * (p + b);
end

function ySolution = naiveODESolver(yDfunc, tspan, y0, zSolution, resolution)
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
        yD = yDfunc(t, y, z);
        if any(isinf(yD)) || any(isnan(yD))
            % assign a random value; who cares
            yD = 1e-2;
        end
        y = y + h * yD;        
        t = t + h;
    end
end