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
points = [1, 1; 0, 0];

sol = initialGenerator(grids, points);
% q1Num = 1 ./ sqrt(1+exp(2*tan(grids)));
% q2Num = 1./sqrt(3*(1./q1Num.^2-1)+1);
% sol = [q1Num; q2Num];

time = 0;
timeStep = 2e-3;
itr = 1;
maxItr = 2e3;
shape = size(sol);

results = zeros(maxItr, 2, numel(grids));
tic
errorV = [];
while itr <= maxItr
    lsol = laplacian(sol, h);
    dsol = derivative(sol, h);
	if adaptive
        sol = sol + (ReactionTerm(sol) + D*(cos(grids).^4.*lsol-2.*sin(grids).*cos(grids).^3.*dsol)) * timeStep;
    else
        % Explicit method
%         sol = sol + (ReactionTerm(sol) + D*lsol) * timeStep;
        % Implicit method
        A = @(v) truncationError(v, shape, timeStep, sol, D, h);
%         sol = reshape(pcg(A, zeros(numel(sol), 1), [], [], [], [], sol), shape(1), shape(2));
        sol = reshape(pcg(A, reshape(sol/timeStep, prod(shape), 1), [], [], [], [], reshape(sol, prod(shape), 1)), shape(1), shape(2));
    end
    
    results(itr, :, :) = sol;
    
    time = time + timeStep;
    itr = itr + 1;
    
    if mod(itr, 50) == 0
        if adaptive
            p = derivative(sol, h) .* cos(grids).^2 - hFunc(sol);
        else
            p = derivative(sol, h) - hFunc(sol);
        end
        hamiltonian = hamiltonianFunc(sol, p);
        fprintf('Iteration % 6d, error in H %.6f.\n', itr, max(abs(hamiltonian)));
        errorV = [errorV max(abs(hamiltonian))];
%         figure(3); hold on; plot(hamiltonian)
    end
end
toc

% figure(2); legend('1', '2', '3', '4', '5')

figure(1)
contour(permute(results(:, 1, :), [1 3 2]))

% figure(3)
% mesh(permute(results(:, 1, :), [1 3 2]))

figure(4); semilogy(errorV)
figure(5); hold on
q1s = 0.05:0.05:0.95;
q2s = 1./sqrt(3*(1./q1s.^2-1)+1);
plot(q1s, q2s, 'r-.')
scatter(sol(1, :), sol(2, :), 'b.')
figure(6); hold on
plot(-2*hFunc(q1s), -2*hFunc(q2s), 'r-.')
scatter(p(1, :), p(2, :), 'b.')

function b = truncationError(v, shape, timeStep, sol, D, h)
    nsol = reshape(v, shape(1), shape(2));
    % Backward Euler
%     t = nsol/timeStep - (ReactionTerm(nsol) + D*laplacian(nsol, h));
    % Crank-Nicolson
    t = nsol/timeStep - ((ReactionTerm(sol) + D*laplacian(sol, h))+(ReactionTerm(nsol) + D*laplacian(nsol, h)))/2;
%     t = (nsol-sol)/timeStep - (ReactionTerm(nsol) + D*laplacian(nsol, h));
    b = reshape(t, numel(nsol), 1);
end

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
%     d(:, 1:end-1) = (v(:, 2:end)-v(:, 1:end-1))/h;
%     d(:, end) = d(:, end-1);
    d(:, 2:end-1) = (v(:, 3:end)-v(:, 1:end-2))/(2*h);
    d(:, end) = 2*d(:, end-1)-d(:, end-2);
    d(:, 1) = 2*d(:, 2)-d(:, 3);
end

function l = laplacian(x, h)
    [m, n] = size(x);
    l = zeros(m, n);
    l(:, 2:end-1) = (x(:, 1:end-2) + x(:, 3:end) - 2*x(:, 2:end-1)) / h^2;
    l(:, 1) = 2*l(:, 2)-l(:, 3);
    l(:, end) = 2*l(:, end-1)-l(:, end-2);
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