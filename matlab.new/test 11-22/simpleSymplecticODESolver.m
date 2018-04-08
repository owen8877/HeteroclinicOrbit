function xSolution = simpleSymplecticODESolver(xDfunc, tspan, x0, ...
    zSolution, resolution, options)
    if ~isfield(options, 'output')
        options.output = true;
    end
    if ~isfield(options, 'method')
        options.method = 'rk';
    end
    
    tStart = tspan(1);
    tEnd = tspan(2);
    lMagic = 1 - 2/sqrt(3);
    hMagic = 1 + 2/sqrt(3);
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
        switch options.method
            case 'euler'
                k = x;
                k = x + h/2 * xDfunc(t, k, z);
                k = x + h/2 * xDfunc(t, k, z);
                x = x + h * xDfunc(t, k, z);
            case 'rk'
                K1 = x + (xDfunc(t, x, z)+xDfunc(t, x, z)*lMagic) * h/4;
                K1 = x + (xDfunc(t, K1, z)+xDfunc(t, x, z)*lMagic) * h/4;
                K1 = x + (xDfunc(t, K1, z)+xDfunc(t, x, z)*lMagic) * h/4;
                K2 = x + (xDfunc(t, K1, z)*hMagic+xDfunc(t, x, z)) * h/4;
                K2 = x + (xDfunc(t, K1, z)*hMagic+xDfunc(t, K2, z)) * h/4;
                K2 = x + (xDfunc(t, K1, z)*hMagic+xDfunc(t, K2, z)) * h/4;
                K1 = x + (xDfunc(t, K1, z)+xDfunc(t, K2, z)*lMagic) * h/4;
                K2 = x + (xDfunc(t, K1, z)*hMagic+xDfunc(t, K2, z)) * h/4;
                
                x = x + h/2 * (xDfunc(t, K1, z) + xDfunc(t, K2, z));
        end
        t = t + h;
        if options.output
            clc
            fprintf('Progress %.1f%%', index / (resolution+1) *100);
        end
    end
end