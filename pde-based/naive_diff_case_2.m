clear; %clc

h0 = 1/16;
adaptive = false;

if adaptive
    h = h0*pi/8;
    grids = (-pi/2):(h0*pi/8):(pi/2);
else
    h = h0;
    grids = (-6+h/2):h:(6-h/2);
end

D = eye(2);
points = [0, 0; 1/2, 0.001];
fFunc = @(x) [4*x.^3 - 1*x; 12*x.^2-1; 24*x];

sol = initialGenerator(grids, points);

time = 0;
exTimeStep = 1e-4; % Explicit
imTimeStep = 1e-3; % Implicit
itr = 1;
maxItr = 2e3;
preExplicit = 100;
shape = size(sol);

for i = 1:3
    figure(i); clf
end

results = zeros(maxItr, 2, numel(grids));
tic
errorV = [];
while itr <= maxItr
%     lsol = laplacian(sol, h);
%     nsol = derivative(sol, h);
%     nh = nhFunc(sol, fFunc);
%     asym = (nh' - nh)*nsol;
%     asym = asymqFunc(sol, nsol, fFunc);
%     sol = sol + (ReactionTerm(sol, fFunc) + D*lsol + asym) * timeStep;
    
    lsol = laplacian(sol, h);
    dsol = derivative(sol, h);
	if adaptive
        sol = sol + (ReactionTerm(sol, fFunc) + D*(cos(grids).^4.*lsol-2.*sin(grids).*cos(grids).^3.*dsol)) * exTimeStep;
    else
        if itr < preExplicit
            % Explicit method
            sol = sol + (ReactionTerm(sol, fFunc) + D*lsol) * exTimeStep;
        else
            % Implicit method
            A = @(v) truncationError(v, shape, imTimeStep, fFunc, sol, D, h);
            sol = reshape(pcg(A, reshape(sol/imTimeStep, prod(shape), 1), [], [], [], [], reshape(sol, prod(shape), 1)), shape(1), shape(2));
        end
    end
    
    results(itr, :, :) = sol;
    
    time = time + imTimeStep;
    itr = itr + 1;
    
    if mod(itr, floor(maxItr/5)) == 0
        figure(2); hold on
        plot(sol(1, :), sol(2, :))
        
        if adaptive
            p = derivative(sol, h) .* cos(grids).^2 - hFunc(sol, fFunc);
        else
            p = derivative(sol, h) - hFunc(sol, fFunc);
        end
        hamiltonV = hamiltonianFunc(sol, p, fFunc);
        figure(3); hold on
        plot(hamiltonV)
    end

    if mod(itr, 100) == 0
        if adaptive
            p = derivative(sol, h) .* cos(grids).^2 - hFunc(sol, fFunc);
        else
            p = derivative(sol, h) - hFunc(sol, fFunc);
        end
        hamiltonV = hamiltonianFunc(sol, p, fFunc);
        fprintf('Iteration % 6d, error in H %.6f.\n', itr, max(abs(hamiltonV)));
        errorV = [errorV max(abs(hamiltonV))];
        
        [m1, n1] = size(sol);
        rhsE = zeros(m1, n1);
        for i = 1:n1
            qV = sol(:, i);
            dqV = dsol(:, i);
            ddqV = lsol(:, i);
            
            nhM = nhFunc(qV, fFunc);
            rhsE(:, i) = ddqV + (nhM' - nhM)*dqV + nlFunc(qV, fFunc) ...
                - nhM*hFunc(qV, fFunc);
        end
    end
end
toc

figure(2); legend('1', '2', '3', '4', '5')
figure(3); legend('1', '2', '3', '4', '5')

figure(1)
contour(permute(results(:, 1, :), [1 3 2]))

% figure(3)
% mesh(permute(results(:, 1, :), [1 3 2]))

function b = truncationError(v, shape, timeStep, fFunc, sol, D, h)
    nsol = reshape(v, shape(1), shape(2));
    % Backward Euler
%     t = nsol/timeStep - (ReactionTerm(nsol, fFunc) + D*laplacian(nsol, h));
    % Crank-Nicolson
    t = nsol/timeStep - ((ReactionTerm(sol, fFunc) + D*laplacian(sol, h))+(ReactionTerm(nsol, fFunc) + D*laplacian(nsol, h)))/2;
    b = reshape(t, numel(nsol), 1);
end

function f = ReactionTerm(q, fFunc)
    f = nlFunc(q, fFunc) - nhhFunc(q, fFunc);
%     [m, n] = size(q);
%     f = zeros(m, n);
%     for i = 1:size(q, 2)
%         qv = q(:, i);
%         f(:, i) = nlFunc(qv, fFunc) - nhFunc(qv, fFunc)' * hFunc(qv, fFunc);
%     end
end

function l = lFunc(q, fFunc)
    fVal = fFunc(q(1, :));
    r = q(2, :)-fVal(1, :);
    l = -2 * r.^2;
end

function nlv = nlFunc(q, fFunc)
    fVal = fFunc(q(1, :));
    r = q(2, :)-fVal(1, :);
    nlv = -4 * r .* [-fVal(2, :); ones(1, size(fVal, 2))];
end

function hv = hFunc(q, fFunc)
    fVal = fFunc(q(1, :));
    hv = 2 * [q(2, :); fVal(2, :).*q(2, :)];
end

function nhthv = nhthFunc(q, fFunc)
    fVal = fFunc(q(1, :));
    nhthv = 4 * [fVal(2, :).*fVal(3, :).*(q(2, :).^2); ...
        q(2, :).*(1+fVal(2, :).^2)];
end

function nhhv = nhhFunc(q, fFunc)
    fVal = fFunc(q(1, :));
    nhhv = 4 * [fVal(2, :).*q(2, :); ...
        q(2, :).*(fVal(3, :).*q(2, :)+fVal(2, :).^2)];
end

function asymv = asymqFunc(q, dq, fFunc)
    fVal = fFunc(q(1, :));
    asymv = 2 * (fVal(3, :).*q(2, :)-1) .* ...
        [dq(2, :); -dq(1, :)];
end

function nhm = nhFunc(q, fFunc)
    fVal = fFunc(q(1, :));
    nhm = 2 * [0 1; fVal(3, :)*q(2, :) fVal(2, :)];
end

function h = hamiltonianFunc(q, p, fFunc)
    h = sum(p.^2/2 + p.*hFunc(q, fFunc) + lFunc(q, fFunc), 1);
%     [m, n] = size(q);
%     h = zeros(1, n);
%     for i = 1:size(q, 2)
%         qv = q(:, i);
%         pv = p(:, i);
%         h(i) = 0.5*(pv'*pv) + pv'*hFunc(qv, fFunc) + lFunc(qv, fFunc);
%     end
end

function d = derivative(v, h)
    [m, n] = size(v);
    d = zeros(m, n);
    d(:, 1:end-1) = (v(:, 2:end)-v(:, 1:end-1))/h;
    d(:, end) = d(:, end-1);
end

function l = laplacian(x, h)
    [m, n] = size(x);
    l = zeros(m, n);
    l(:, 2:end-1) = (x(:, 1:end-2) + x(:, 3:end) - 2*x(:, 2:end-1)) / h^2;
    l(:, 1) = l(:, 2);
    l(:, end) = l(:, end-1);
end

function sol = initialGenerator(grid, points)
    d1 = round(numel(grid)/6); d2 = round(numel(grid)*5/6);
    grid_trunc = grid(d1+1:d2); grid_trunc = (grid_trunc - grid_trunc(1)) / (grid_trunc(end) - grid_trunc(1));
    state_gap = (points(2, :) - points(1, :))';
    sol = [ ...
        points(1, :)' .* ones(1, d1), ...
        points(1, :)' + state_gap .* ((grid_trunc.^3) .* (10-15*grid_trunc+6*grid_trunc.^2)), ...
        points(2, :)' .* ones(1, numel(grid)-d2), ...
    ];
end