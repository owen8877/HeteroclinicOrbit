function xSolution = naiveODESolver(xDfunc, tspan, x0, zSolution, resolution, charFunc)
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
        if index == 1
            xD = charFunc(x, z);
        else
            xD = xDfunc(t, x, z);
        end
        x = x + h * xD;        
        t = t + h;
    end
end