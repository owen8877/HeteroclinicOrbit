function ySolution = myODESolver(yDfunc, tspan, y0, zSolution, round, varargin)
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
        
        hint = 1;
        
        if useZ
%             z = interp1(tArray_zSolution, zArray_zSolution, t, 'linear', 'extrap')';
%             z_ = interp1(tArray_zSolution, zArray_zSolution, t+h, 'linear', 'extrap')';
            [z, hint] = quickInterpolate(tArray_zSolution, zArray_zSolution, t, hint, round);
            z = z';
            [z_, ~] = quickInterpolate(tArray_zSolution, zArray_zSolution, t+h, hint, round);
            z_ = z_';
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
