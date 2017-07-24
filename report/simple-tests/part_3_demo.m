clear;
clc;

B = [[-1 4]; [-4 -1]];
p0 = [4; 1];
x1 = [0; 0];

drawInterval = 20;
iterationLimit = 100;
preciseOn = true;

t = linspace(0, 10, 1001);
pPrecise = zeros(2, size(t, 2));
xPrecise = zeros(2, size(t, 2));
for i = 1:size(t, 2)
    pPrecise(:, i) = expm(B' * t(i)) * p0;
    xPrecise(:, i) = (B + B') \ pPrecise(:, i);
end

for depth = 4:4
    iteration = 0;
    p1Collection = [];
    previousP1 = [];

    resolution = 10^depth;
    criterion = min(1e-6, 1e-2/resolution);
    sArray = 0:1/resolution:1;
    % init xSolution
    xSolution = [1-sArray; sArray];
    C = sqrt(2);

    %% Iteration
    while true
        pSolution = naiveODESolver(@(s, p, x) pDfunc(s, p, x, B, C), [0 1], p0, xSolution, resolution);
        delta = 1 / (iteration+1);
        xSolution = naiveODESolver(@(s, x, p) xDfunc(s, x, p, B, C), [1 0], x1 + [delta; 0], pSolution, resolution);

        iteration = iteration + 1;
        if iteration - fix(iteration / drawInterval) * drawInterval == 0
            figure;
            s1 = subplot(1,2,1); hold(s1, 'on');
            plot(pSolution(1, :), pSolution(2, :), 'b-', 'DisplayName', 'numerical');
            plot(pPrecise(1, :), pPrecise(2, :), 'r-', 'DisplayName', 'analytical');
            
            s2 = subplot(1,2,2); hold(s2, 'on');
            plot(xSolution(1, :), xSolution(2, :), 'b-', 'DisplayName', 'numerical');
            plot(xPrecise(1, :), xPrecise(2, :), 'r-', 'DisplayName', 'analytical');
            legend('show');
            title(['Iteration ' num2str(iteration)])
        end

        %% update C
        C = 0;
        for i = 1:(size(xSolution, 2) - 1)
            C = C + norm(xSolution(:, i) - xSolution(:, i+1));
        end
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

    %save(['iteration-' num2str(depth) '.mat'], 'xSolution', 'pSolution', 'p1Collection');
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