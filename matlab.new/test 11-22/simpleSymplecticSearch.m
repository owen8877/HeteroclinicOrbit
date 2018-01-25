function xSolution = simpleSymplecticSearch(xDfunc, origin, target, step, mode, indicator)
    h = step;
    t = 0; index = 1;
    
    x = origin;
    xSolution = zeros(size(origin, 1), max(2/step, 1e4));
    nearFlag = false;
    modifiedLastDistance = -Inf;
    mu = 1;
    if mode == 0 % direction
        indicator = indicator / norm(indicator);
    end
    
    while true
        xSolution(:, index) = x;
        
        k = x;
        k = x + h/2 * xDfunc(t, k);
        k = x + h/2 * xDfunc(t, k);
        x = x + h * xDfunc(t, k);
        t = t + h;
        index = index + 1;
        
        % see close enough
        if norm(target-x) < abs(h) * 1000
            if mode == 0 % direction
                tmpDistance = norm(target-x) + mu * res(x, indicator);
            elseif mode == 1 % decomp matrix
                tmpDistance = indicator' * (target-x);
                tmpDistance = abs(tmpDistance) / norm(tmpDistance);
                tmpDistance = tmpDistance(1);
            end
            
            if ~nearFlag
                nearFlag = true;
                modifiedLastDistance = tmpDistance;
            elseif modifiedLastDistance > tmpDistance
            % elseif tmpDistance < 0.3
                xSolution = xSolution(:, 1:index-1);
                break
            else
                modifiedLastDistance = tmpDistance;
            end
        end
                
%         clc
%         fprintf('Distance %.4e', norm(target-x));
        
        if norm(target-x) > 3
            xSolution = xSolution(:, 1:index-1);
            break
        end
    end
    fprintf('\n')
end

function r = res(x, rdelta)
    p = x' * rdelta;
    r = sqrt(x'*x - p^2);
end