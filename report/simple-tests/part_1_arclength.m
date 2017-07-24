clear;
clc;

B = [[-1 4]; [-4 -1]];
p0 = [4; 1];
x1 = [0; 0];

drawInterval = nan;
iterationLimit = nan;

xError = [];
depthRange = 2:6;
for depth = depthRange
    iteration = 0;
    p1Collection = [];
    previousP1 = [];

    resolution = 10^depth;
    criterion = 1e-7;
    sArray = 0:1/resolution:1;
    tArray = sArray * 0;
    % init xSolution
    xSolution = [1-sArray; sArray];
    C = sqrt(2);

    %% Iteration
    while true
        pSolution = naiveODESolver(@(s, p, x) pDfunc(s, p, x, B, C), [0 1], p0, xSolution, resolution);
        delta = 9.993063e-01;
        xSolution = naiveODESolver(@(s, x, p) xDfunc(s, x, p, B, C), [1 0], x1 + [delta; 0], pSolution, resolution);

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
        
%         %% update tArray
%         for i = 1:resolution
%             b = - B*xSolution(:, i);
%             lambdar = C / norm(pSolution(:, i)+b);
%             tArray(i+1) = tArray(i) + lambdar / resolution;
%         end

        %% xError is just simply the distance between the 'ending' points
        %  of xSolution -- in fact, it is xSolution(0).
        xError(depth, iteration) = norm(xSolution(:, 1) - (B + B') \ p0);

        %% output iteration info
        clc
        fprintf('(Depth %d) Iteration %d completed.\n', depth, iteration);

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
    
    % save(['iteration-' num2str(depth) '.mat'], 'xSolution', 'pSolution', 'p1Collection');
end

figure
hold('on');
for depth = depthRange
    plot(1:size(xError, 2), log10(xError(depth, :)), 'DisplayName', sprintf('1e%d', depth));
end
title(sprintf('log10(xError)\ncriterion %.1e', criterion));
legend('show');

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