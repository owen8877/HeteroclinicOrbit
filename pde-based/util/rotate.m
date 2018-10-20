function newq = rotate(q)
    [n, m] = size(q); clusterN = n / 2;
    newq = zeros(n, m);
    newq(:, 1) = q(:, 1);
    
    sumtheta = 0;
    for i = 2:m
        new = reshape(q(:, i), 2, clusterN);
        old = reshape(q(:, i-1), 2, clusterN);
        I = norm(q(:, i))^2;
        displacement = new - old;
        angleDisplacement = sum(cross(old, displacement));
        
        sumtheta = sumtheta + angleDisplacement / I;
        
        newq(:, i) = reshape(orthRotate(new, -sumtheta, clusterN), n, 1);
    end
end

function v = orthRotate(v, theta, clusterN)
    for i = 1:clusterN
        v(:, i) = [cos(theta) -sin(theta); sin(theta) cos(theta)] * v(:, i);
    end
end

function c = cross(A, B)
   c = A(1, :) .* B(2, :) - A(2, :) .* B(1, :);
end