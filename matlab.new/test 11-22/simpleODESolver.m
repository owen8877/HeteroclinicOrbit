function xSolution = simpleODESolver(xDfunc, tspan, x0, zSolution, resolution)
    tStart = tspan(1);
    tEnd = tspan(2);
    h = (tEnd - tStart) / resolution;
    t = tStart;
    directionPos = tEnd > tStart;
    
    x = x0;
    xSolution = zeros(size(x0, 1), resolution+1);
    
    for index = 1:(resolution+1)
        if directionPos
            realIndex = index;
        else
            realIndex = resolution + 2 - index;
        end
        
        % record x
        xSolution(:, realIndex) = x;
        
        z = zSolution(:, realIndex);
        xD = xDfunc(t, x, z);
        x = x + h * xD;        
        t = t + h;
    end
end