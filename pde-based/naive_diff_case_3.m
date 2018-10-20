clear; %clc

h = 1/64;
grids = (-4+h/2):h:(4-h/2);

global weight gamma
dimension = 2;
D = eye(dimension);
weight = 1:dimension;
gamma = @(x) (x+1).*(x-1);
points = [-1; 1] .* ones(1, dimension);

sol = initialGenerator(grids, points);

time = 0;
timeStep = 1e-4;
itr = 1;
maxItr = 10000;

results = zeros(maxItr*2, dimension, numel(grids));
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
%         figure(5); hold on; plot(hamiltonian);
%         figure(6); hold on; plot3(sol(1, :), sol(2, :), sol(3, :)); grid on
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
figure(5); plot(hamiltonian);
% figure(6); hold on; scatter3(sol(1, :), sol(2, :), sol(3, :), 'b.'); grid on
figure(6); hold on; scatter(sol(1, :), sol(2, :), 'b.'); grid on

function f = ReactionTerm(q)
    f = nlFunc(q) - nhhFunc(q);
end

function l = lFunc(q)
    l = zeros(1, size(q, 2));
end

function nlv = nlFunc(q)
    nlv = q * 0;
end

function nhhv = nhhFunc(q)
    global weight gamma
    G = gamma(q);
    Gp1 = circshift(G, [1, 0]);
    Gn1 = circshift(G, [-1, 0]);
    weightp1 = circshift(weight, [0, 1]);
    qp1 = circshift(q, [1, 0]);
    qn1 = circshift(q, [-1, 0]);
    
    t = weight' .* Gn1 .^ 2 + weightp1' .* Gp1 .^ 2;
    tp1 = circshift(t, [1, 0]);
    tn1 = circshift(t, [-1, 0]);
    
    nhhv = 64 * (qp1.^2 .* q .* Gp1.^2 .* G .* weightp1' .* tp1) ...
        + 16 * q .* G .* (G + 2 .* q.^2) .* t.^2 ...
        + 64 * (qn1.^2 .* q .* Gn1.^2 .* G .* weight' .* tn1);
end

function hv = hFunc(q)
    global weight gamma
    G = gamma(q);
    Gp1 = circshift(G, [1, 0]);
    Gn1 = circshift(G, [-1, 0]);
    weightp1 = circshift(weight, [0, 1]);
    hv = 4 * q .* G .* (weight'.* Gn1.^2 + weightp1' .* Gp1.^2);
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