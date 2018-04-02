function solution = pulling(oldSolution, threshold, beta, target)
    % pulling starts at (1-threshold), with param beta, to the target point
    n = size(oldSolution, 2);
    solution = oldSolution;
    
    leftPoint = floor((1-threshold)*n);
    grids = ((leftPoint:n)-leftPoint)/(n-leftPoint);
    lambda = grids .^ beta;
    solution(:, leftPoint:n) = solution(:, leftPoint:n) + lambda .* (target - solution(:, leftPoint:n));
end

