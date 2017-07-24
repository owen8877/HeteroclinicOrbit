clear;
clc;

load('high-order-equation-B.mat');
% B = [[-1 4]; [-4 -1]];
p0 = [1; 2; 3];
x1 = [0; 0; 0];

drawInterval = nan;
iterationLimit = nan;
preciseOn = true;
resolutionFactor = 10;

t = linspace(0, 20, 1001);
pPrecise = zeros(3, size(t, 2));
xPrecise = zeros(3, size(t, 2));
for i = 1:size(t, 2)
    pPrecise(:, i) = expm(B' * t(i)) * p0;
    xPrecise(:, i) = (B + B') \ pPrecise(:, i);
end

for depth = 5:5
    iteration = 0;
    p1Collection = [];
    previousP1 = [];
    xInitInfo = [];

    resolution = 10^depth;
    criterion = 1e-1 / resolution;
    sArray = 0:1/resolution:1;
    % init xSolution
    xSolution = [1-sArray; sArray; 0 * sArray];
    C = sqrt(2);
    xInitial = 1/4;

    %% Iteration
    while true
        pSolution = naiveODESolver(@(s, p, x) pDfunc(s, p, x, B, C), [0 1], p0, xSolution, resolution);
        xSolution = naiveODESolver(@(s, x, p) xDfunc(s, x, p, B, C), [1 0], x1 + [xInitial; 0; 0], pSolution, resolution);

        iteration = iteration + 1;
        
        if iteration - fix(iteration / drawInterval) * drawInterval == 0
            figure;
            s1 = subplot(1,2,1); hold(s1, 'on');
            plot3(pSolution(1, :), pSolution(2, :), pSolution(3, :), 'b-', 'DisplayName', 'numerical');
            plot3(pPrecise(1, :), pPrecise(2, :), pPrecise(3, :), 'r-', 'DisplayName', 'analytical');
            legend('show');
            title(['pSolution Iteration ' num2str(iteration)])
            
            s2 = subplot(1,2,2); hold(s2, 'on');
            plot3(xSolution(1, :), xSolution(2, :), xSolution(3, :), 'b-', 'DisplayName', 'numerical');
            plot3(xPrecise(1, :), xPrecise(2, :), xPrecise(3, :), 'r-', 'DisplayName', 'analytical');
            legend('show');
            title(['xSolution Iteration ' num2str(iteration)])
        end

        %% update C
        C = 0;
        for i = 1:(size(xSolution, 2) - 1)
            C = C + norm(xSolution(:, i) - xSolution(:, i+1));
        end
        clc
        fprintf('Parameters:\n');
        fprintf('\tResolution\t%.0e\n', resolution);
        fprintf('\txInitial\t%.4e\n', xInitial);
        fprintf('\tcriterion\t%.0e\n', criterion);
        fprintf('Iteration %d completed.\n', iteration);

        if iteration > iterationLimit
            break
        end

        %% watch if p(1) converges
        if iteration ~= 1
            p1 = pSolution(:, end);
            endDistance = norm(p1 - previousP1);
            if endDistance < criterion
                xInitial = xInitial / 2;
                if xInitial < resolutionFactor/resolution
                    break
                end
            end
            fprintf('probable error in p: %e\n', endDistance);
        end
        xInitInfo(iteration) = xInitial;
        previousP1 = pSolution(:, end);
        p1Collection(:, iteration) = previousP1;
    end
    
    tArray = sArray * 0;
    %% update tArray
    for i = 1:resolution
        p = pSolution(:, i);
        b = - B*xSolution(:, i);
        tArray(i+1) = tArray(i) + C / norm(p+b) * 1 / resolution;
    end

    pP = zeros(3, size(tArray, 2));
    xP = zeros(3, size(tArray, 2));
    for i = 1:size(tArray, 2)
        pP(:, i) = expm(B' * tArray(i)) * p0;
        xP(:, i) = (B + B') \ pP(:, i);
    end

    pDiff = pSolution - pP;
    xDiff = xSolution - xP;

    pErr = zeros(1, size(tArray, 2));
    xErr = zeros(1, size(tArray, 2));
    for i = 1:size(tArray, 2)
        pErr(i) = norm(pDiff(:, i));
        xErr(i) = norm(xDiff(:, i));
    end
    
    fprintf('actual error in p: %e\n', max(pErr));

    save(['s-iteration-' num2str(depth) '.mat'], 'xSolution', 'pSolution', 'p1Collection', 'xInitInfo', 'pErr', 'xErr');
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