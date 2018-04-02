function xSolution = simpleSymplecticODESolver(xDfunc, tspan, x0, zSolution, resolution)
    % xDfunc(t, x, z), tspan, x0, zSolution, resolution
    tStart = tspan(1);
    tEnd = tspan(2);
    h = (tEnd - tStart) / resolution;
    t = tStart;
    directionPos = tEnd > tStart;
    
    x = x0;
    xSolution = zeros(size(x0, 1), round(resolution+1));
    
    ignoreZ = size(zSolution, 2) < resolution;
    
    for index = 1:(resolution+1)
        if directionPos
            realIndex = index;
        else
            realIndex = round(resolution + 2) - index;
        end
        
        % record x
        xSolution(:, realIndex) = x;
        
        if ignoreZ
            z = 0;
        else
            z = zSolution(:, realIndex);
        end
        k = x;
        k = x + h/2 * xDfunc(t, k, z);
        k = x + h/2 * xDfunc(t, k, z);
        x = x + h * xDfunc(t, k, z);
        t = t + h;
        clc
        fprintf('Progress %.1f%%', index / (resolution+1) *100);
    end
end