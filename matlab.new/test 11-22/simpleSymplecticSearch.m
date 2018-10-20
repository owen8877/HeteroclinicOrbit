function xSolution = simpleSymplecticSearch(xDfunc, origin, target, step, ...
        mode, indicator, options)
    if ~isfield(options, 'output')
        options.output = true;
    end
    if ~isfield(options, 'endCriterion')
        options.endCriterion = 40;
    end
    if ~isfield(options, 'method')
        options.method = 'euler';
    end
    if ~isfield(options, 'torlerance')
        options.torlerance = 1e-15;
    end
    
    h = step;
    t = 0; index = 1;
    
    x = origin;
    xSolution = zeros(size(origin, 1), max(floor(abs(2/step)), 1e4));
    nearFlag = -1; % -1 from beginning, 0 for outside equilibrium, 1 for near
    modifiedLastDistance = -Inf;
    mu = 10;
    if mode == 0 % direction
        indicator = indicator / norm(indicator);
    end
    
    lMagic = 1 - 2/sqrt(3);
    hMagic = 1 + 2/sqrt(3);
    
    while true
        xSolution(:, index) = x;
        
        switch options.method
            case 'euler'
                newk = newtonSubSolver(x, @(k) x+h/2*xDfunc(t, k));
                x = x + h * xDfunc(t, newk);
            case 'rk'
%                 K1 = x + (xDfunc(t, x)+xDfunc(t, x)*lMagic) * h/4;
%                 K2 = x + (xDfunc(t, K1)*hMagic+xDfunc(t, x)) * h/4;
%                 K1 = x + (xDfunc(t, K1)+xDfunc(t, K2)*lMagic) * h/4;
%                 K2 = x + (xDfunc(t, K1)*hMagic+xDfunc(t, K2)) * h/4;
%                 K1 = x + (xDfunc(t, K1)+xDfunc(t, K2)*lMagic) * h/4;
%                 K2 = x + (xDfunc(t, K1)*hMagic+xDfunc(t, K2)) * h/4;

                K1 = x + (xDfunc(t, x)+xDfunc(t, x)*lMagic) * h/4;
                K1 = x + (xDfunc(t, K1)+xDfunc(t, x)*lMagic) * h/4;
                K1 = x + (xDfunc(t, K1)+xDfunc(t, x)*lMagic) * h/4;
                K2 = x + (xDfunc(t, K1)*hMagic+xDfunc(t, x)) * h/4;
                K2 = x + (xDfunc(t, K1)*hMagic+xDfunc(t, K2)) * h/4;
                K2 = x + (xDfunc(t, K1)*hMagic+xDfunc(t, K2)) * h/4;
                K1 = x + (xDfunc(t, K1)+xDfunc(t, K2)*lMagic) * h/4;
                K2 = x + (xDfunc(t, K1)*hMagic+xDfunc(t, K2)) * h/4;
                
                x = x + h/2 * (xDfunc(t, K1) + xDfunc(t, K2));
        end
        
        t = t + h;
        index = index + 1;
        
        % see close enough
%         if norm(target-x) < abs(h) * 100
        if norm(target-x) < 0.1
            if mode == 0 % direction
                tmpDistance = norm(target-x) + mu * res(x, indicator);
            elseif mode == 1 % decomp matrix
                tmpDistance = indicator' * (target-x);
                tmpDistance = abs(tmpDistance) / norm(tmpDistance);
                tmpDistance = -tmpDistance(1);
            elseif mode == 2 % distance mode
                tmpDistance = norm(target-x);
            end
            
            if nearFlag == 0
                nearFlag = 1;
                modifiedLastDistance = tmpDistance;
            elseif nearFlag == 1
                if modifiedLastDistance < tmpDistance
                % elseif tmpDistance < 0.3
                    xSolution = xSolution(:, 1:index-1);
                    break
                else
                    modifiedLastDistance = tmpDistance;
                end
            end
        else
            if nearFlag == -1
                nearFlag = 0;
            end
        end
        
        if options.output && mod(index, 1000) == 0
            clc
            fprintf('Distance %.4e\t Index % 6d', norm(target-x), index);
        end
        
        if norm(target-x) > options.endCriterion || abs(index * step) > 20
            xSolution = xSolution(:, 1:index-1);
            break
        end
    end
    
    if options.output
        fprintf('\n')
    end
end

function r = res(x, rdelta)
    p = x' * rdelta;
    r = sqrt(x'*x - p^2);
end