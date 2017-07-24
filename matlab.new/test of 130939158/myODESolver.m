function ySolution = myODESolver(yDfunc, tspan, y0, zSolution, resolution, varargin)
    useZ = true;
    if nargin < 4
        zSolution = [];
        resolution = 1e3;
        useZ = false;
    end
    
    tStart = tspan(1);
    tEnd = tspan(2);
    h = (tEnd - tStart) / resolution;
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
    ySolution.t = zeros(1, resolution+1);
    ySolution.y = zeros(size(y0, 1), resolution+1);
    
    hint = 1;
    
    while true
        % record t and y
        ySolution.t(index) = t;
        ySolution.y(:, index) = y;
        
        if useZ
%             z = interp1(tArray_zSolution, zArray_zSolution, t, 'linear', 'extrap')';
%             z_ = interp1(tArray_zSolution, zArray_zSolution, t+h, 'linear', 'extrap')';
            [zt, hint] = quickInterpolate(tArray_zSolution, zArray_zSolution, t, hint);
            z = zt';
%             [z_, ~] = quickInterpolate(tArray_zSolution, zArray_zSolution, t+h, hint);
%             z_ = z_';
%             y_ = y + h * yDfunc(t, y, z);
%             y = y + h/2 * (yDfunc(t, y, z) + yDfunc(t+h, y_, z_));
            yD = yDfunc(t, y, z);
            if any(isinf(yD)) || any(isnan(yD))
                % assign a random value; who cares
                yD = 1e-2;
            end
%             y_ = y + h * yD;
%             yD = yDfunc(t, y_, z);
%             if any(isinf(yD)) || any(isnan(yD))
%                 % assign a random value; who cares
%                 yD = 1e-2;
%             end
            y = y + h * yD;
        else
%             y_ = y + h * yDfunc(t, y);
%             y = y + h/2 * (yDfunc(t, y) + yDfunc(t+h, y_));
            yD = yDfunc(t, y);
            if any(isinf(yD)) || any(isnan(yD))
                % assign a random value; who cares
                yD = 1e-2;
            end
%             y_ = y + h * yD;
%             yD = yDfunc(t, y_);
%             if any(isinf(yD)) || any(isnan(yD))
%                 % assign a random value; who cares
%                 yD = 1e-2;
%             end
            y = y + h * yD;
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
