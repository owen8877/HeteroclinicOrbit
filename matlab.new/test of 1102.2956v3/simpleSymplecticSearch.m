function xSolution = simpleSymplecticSearch(xDfunc, origin, target, step, probableSize, rdelta)
    h = step;
    t = 0; index = 1;
    
    x = origin;
    xSolution = zeros(size(origin, 1), max(probableSize, 1e4));
    nearFlag = false;
    modifiedLastDistance = 0;
    mu = 10;

    while true
        xSolution(:, index) = x;
        
        k = x;
        k = x + h/2 * xDfunc(t, k);
        k = x + h/2 * xDfunc(t, k);
        x = x + h * xDfunc(t, k);
        t = t + h;
        index = index + 1;
        
        % see close enough
        if norm(target-x) < h * 100
            if ~nearFlag
                nearFlag = true;
                modifiedLastDistance = norm(target-x) + mu * res(x, rdelta);
            elseif modifiedLastDistance < norm(target-x) + mu * res(x, rdelta)
                xSolution = xSolution(:, 1:index-1);
                break
            else
                modifiedLastDistance = norm(target-x) + mu * res(x, rdelta);
            end
        end
        
        clc
        fprintf('Distance %.4e', norm(target-x));
    end
    fprintf('\n')
end

function r = res(x, rdelta)
    p = x' * rdelta;
    r = sqrt(x'*x - p^2);
end