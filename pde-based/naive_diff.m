clear; %clc

h = 1/32;
grid = (-4+h/2):h:(4-h/2);
D = eye(2);
points = [0, 0; 1, 1];

sol = initialGenerator(grid, points);

time = 0;
timeStep = 1e-4;
itr = 1;
maxItr = 5e3;

results = zeros(maxItr*2, 2, numel(grid));
tic
errorV = [];
while itr <= maxItr
    lsol = laplacian(sol, h);
    sol = sol + (ReactionTerm(sol) + D*lsol) * timeStep;
    
    results(itr, :, :) = sol;
    
    time = time + timeStep;
    itr = itr + 1;
    
    if mod(itr, 100) == 0
        p = derivative(sol, h) - hFunc(sol);
        hamiltonian = hamiltonianFunc(sol, p);
        fprintf('Iteration % 6d, error in H %.6f.\n', itr, max(abs(hamiltonian)));
        errorV = [errorV max(abs(hamiltonian))];
    end
end
toc

tic
itr = 0;
while itr <= maxItr
    lsol = laplacian(sol, h);
    sol = sol + (ReactionTerm(sol) + D*lsol) * timeStep * 2;
    
    results(itr+maxItr, :, :) = sol;
    
    time = time + timeStep;
    itr = itr + 1;
    
    if mod(itr, 1) == 0
        p = derivative(sol, h) - hFunc(sol);
        hamiltonian = hamiltonianFunc(sol, p);
        fprintf('Iteration % 6d, error in H %.6f.\n', itr, max(abs(hamiltonian)));
        errorV = [errorV max(abs(hamiltonian))];
    end
end
toc

% figure(2); legend('1', '2', '3', '4', '5')
% 
% figure(1)
% contour(permute(results(:, 1, :), [1 3 2]))

% figure(3)
% mesh(permute(results(:, 1, :), [1 3 2]))

figure(4); plot(errorV)

% p = derivative(sol, h) - H(sol);
% hamiltonV = hamilton(sol, p);
% figure(4); plot(hamiltonV);

% range = 110:180;
% aa = a(range);
% bb = b(range);
% gridRes = 32; gridMax = 16;
% 
% xq = linspace(1/gridRes, gridMax, gridRes * gridMax);
% aa = interp1((range-range(1))*h, aa, xq, 'linear', 'extrap');
% bb = interp1((range-range(1))*h, bb, xq, 'linear', 'extrap');
% 
% aaa = [flip(aa), aa];
% bbb = [flip(bb), bb];

% save(sprintf('init-%.2f/init_a.dat', k(9)), 'aaa', '-ascii');
% save(sprintf('init-%.2f/init_b.dat', k(9)), 'bbb', '-ascii');

% figure
% subplot(1, 2, 1)
% contour(resultsa); colorbar
% subplot(1, 2, 2)
% contour(resultsb); colorbar

function f = ReactionTerm(q)
    f = nlFunc(q) - nhhFunc(q);
end

function l = lFunc(q)
    L = 3;
    q1 = q(1, :); q2 = q(2, :);
    r = q1.^2.*(1-q2.^2) - L*q2.^2.*(1-q1.^2);
    l = r.^2;
end

function nlv = nlFunc(q)
    L = 3;
    q1 = q(1, :); q2 = q(2, :);
    l = lFunc(q);
    nlv = 2*l .* [2*q1.*(1-q2.^2)+2*L.*q1.*q2.^2; -2*L.*q2.*(1-q1.^2)-2*q2.*q1.^2];
end

function nhhv = nhhFunc(q)
    nhhv = (q-q.^3) .* (1-3*q.^2);
end

function hv = hFunc(q)
    hv = q-q.^3;
end

function nhm = nhFunc(q)
    nhm = diag(1-3*q.^2);
end

function h = hamiltonianFunc(q, p)
    h = 0.5 * sum(p.^2, 1) + sum(p.*hFunc(q), 1) + lFunc(q);
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