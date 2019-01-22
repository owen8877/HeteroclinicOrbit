clear; %clc
addpath util

clusterN = 7;

s3d2 = sqrt(3)/2; s62 = 2^(1/6);
positionA = [ ...
    0, 0; ...
    s62, 0; ...
    s62/2, s62*s3d2; ...
    -s62/2, s62*s3d2; ...
    -s62, 0; ...
    -s62/2, -s62*s3d2; ...
    s62/2, -s62*s3d2; ...
]'; % State A
positionB = [ ...
    0, 0; ...
    s62*1.5, -s62*s3d2; ...
    s62, 0; ...
    s62/2, s62*s3d2; ...
    -s62, 0; ...
    -s62/2, -s62*s3d2; ...
    s62/2, -s62*s3d2; ...
]'; % State B
positionA = centerize(positionA);
positionB = centerize(positionB);

options.maxitr = 100;
options.beta = 0.5;
options.alpha = 0.5;
positionA = initialSolver(positionA, options);
positionB = initialSolver(positionB, options);

%% React-Diff Equation
h0 = 1/8;
adaptive = false;

if adaptive
    h = h0*pi/8;
    grids = (-pi/2):(h0*pi/8):(pi/2);
else
    h = h0;
    grids = (-4+h/2):h:(4-h/2);
end

points = [reshape(positionB, 1, 2*clusterN); reshape(positionA, 1, 2*clusterN)];

if exist('backup/matlab.mat', 'file')
    load('backup/matlab.mat');
else
    sol = initialGenerator(grids, points);
end

time = 0;
exTimeStep = 1e-5;
imTimeStep = 1e-4;
itr = 1;
maxItr = 2000;
preExp = 2000;
shape = size(sol);

for i = 1:4
    figure(i); clf
end

results = zeros(maxItr, 2*clusterN, numel(grids));
tic
errorV = [];
times = [];
while itr <= maxItr
    lsol = laplacian(sol, h);
    dsol = derivative(sol, h);
	if adaptive
        timeStep = exTimeStep;
        sol = sol + (ReactionTerm(sol) + (cos(grids).^4.*lsol-2.*sin(grids).*cos(grids).^3.*dsol)) * timeStep;
    else
        if itr < preExp
            timeStep = exTimeStep;
            % Explicit method
            sol = sol + (ReactionTerm(sol) + lsol) * exTimeStep;
        else
            timeStep = imTimeStep;
            % Implicit method
            guessSol = sol + (ReactionTerm(sol) + lsol) * exTimeStep;
            A = @(v) truncationError(v, shape, timeStep, sol, h);
            sol = reshape(pcg(A, reshape(sol/timeStep, prod(shape), 1), [], [], [], [], reshape(guessSol, prod(shape), 1)), shape(1), shape(2));
        end
    end
    
    results(itr, :, :) = sol;
    
    time = time + timeStep;
    itr = itr + 1;
    
    if itr < 20 || mod(itr, 20) == 0
        for i = 1:size(sol, 2)
            hV(:, i) = hFunc(sol(:, i));
        end
        if adaptive
            p = derivative(sol, h) .* cos(grids).^2 - hV;
        else
            p = derivative(sol, h) - hV;
        end
        hamiltonV = hamiltonianFunc(sol, p);
        fprintf('Iteration: %d, |H|: %.6f\n', itr, max(abs(hamiltonV)));
        errorV = [errorV max(abs(hamiltonV))];
        times = [times time];
    end
    
    if mod(itr, floor(maxItr/5)) == 0 && false
        figure(2); hold on
        potentials = zeros(size(sol, 2), 1);
        arclength = zeros(size(sol, 2), 1);
        for i = 1:size(sol, 2)
            potentials(i) = potential(reshape(sol(:, i), 2, clusterN));
            hV(:, i) = hFunc(sol(:, i));
            if i == size(sol, 2)
                continue
            end
            arclength(i) = norm(sol(:, i+1)-sol(:, i));
        end
        potentials(end) = potentials(end-1);
        plot(cumsum(arclength), potentials*4)
        
        figure(3); hold on
        plot(hamiltonV)
        
        figure(4); hold on
        potentialPreps = zeros(size(sol, 2), 1);
        for i = 1:size(sol, 2)-1
            nh = nhFunc(sol(:, i));
            aaa = -nh'*p(:, i);
            d = p(:, i+1)-p(:, i);
            d = d/norm(d);
            potentialPreps(i) = (aaa'*d)/norm(aaa);
        end
        plot(cumsum(arclength), potentialPreps)
    end

%     if mod(itr, 100) == 0
%         p = derivative(sol, h) - hFunc(sol);
%         hamiltonV = hamiltonianFunc(sol, p);
%         fprintf('Iteration % 6d, error in H %.6f.\n', itr, max(abs(hamiltonV)));
%         errorV = [errorV max(abs(hamiltonV))];
%         
%         [m1, n1] = size(sol);
%         rhsE = zeros(m1, n1);
%         for i = 1:n1
%             qV = sol(:, i);
%             dqV = nsol(:, i);
%             ddqV = lsol(:, i);
%             
%             nhM = nhFunc(qV);
%             rhsE(:, i) = ddqV - nhM*hFunc(qV);
%         end
%     end
end
toc

% save('backup/matlab.mat', 'sol')

% figure(2); hold on; legend('1', '2', '3', '4', '5')
% figure(3); hold on; legend('1', '2', '3', '4', '5')

figure(1)
contour(permute(results(:, 1, :), [1 3 2]))

function b = truncationError(v, shape, timeStep, sol, h)
    nsol = reshape(v, shape(1), shape(2));
    % Backward Euler
    t = nsol/timeStep - (ReactionTerm(nsol) + laplacian(nsol, h));
    % Crank-Nicolson
%     t = nsol/timeStep - ((ReactionTerm(sol) + laplacian(sol, h))+(ReactionTerm(nsol) + laplacian(nsol, h)))/2;
%     t = (nsol-sol)/timeStep - (ReactionTerm(nsol) + laplacian(nsol, h));
    b = reshape(t, numel(nsol), 1);
end

function f = ReactionTerm(q)
%     f = nlFunc(q, fFunc) - nhhFunc(q, fFunc);
    [m, n] = size(q);
    f = zeros(m, n);
    for i = 1:n
        qv = q(:, i);
        f(:, i) = - nhFunc(qv)' * hFunc(qv);
    end
end

function hv = hFunc(q)
    tc = numel(q);
    position = reshape(q, 2, tc/2);
    hv = reshape(g_potential(position), tc, 1);
end

function nhm = nhFunc(q)
    h = 1e-8;
    tc = numel(q);
    
    nhm = zeros(tc);
    for i = 1:tc
        e = zeros(tc, 1); e(i) = h;
        nhm(:, i) = (hFunc(e+q)-hFunc(q))/h;
    end
end

function h = hamiltonianFunc(q, p)
%     h = sum(p.^2/2 + p.*hFunc(q));
    [~, n] = size(q);
    h = zeros(1, n);
    for i = 1:n
        qv = q(:, i);
        pv = p(:, i);
        h(:, i) = pv'*pv/2 + pv'*hFunc(qv);
    end
end

%% Initial Boundary
function position = centerize(position)
    center = mean(position, 2);
    position = position - center;
end

function V = potential(position)
    n = size(position, 2);
    d = distance(position);
    
    halfd = reshape(triu(d, 1), n^2, 1);
    halfd = halfd(halfd ~= 0);
    
    V = sum(halfd.^-12-halfd.^-6);
end

function g = g_potential(position)
    n = size(position, 2);
    d = distance(position);
    
    g = zeros(2, n);
    for l = 1:2
        for i = 1:n
            tmp = 0;
            for j = 1:n
                if i == j
                    continue
                end
                tmp = tmp + (-12*d(i, j)^-13+6*d(i, j)^-7)*(position(l, i)-position(l, j))/d(i, j);
            end
            g(l, i) = tmp;
        end
    end
end

function d = distance(position)
    n = size(position, 2);
    d = zeros(n);
    for i = 1:n
        for j = i+1:n
            d(i, j) = norm(position(:, i)-position(:, j));
        end
    end
    d = d + d';
end

function position = initialSolver(position, options)
    for itr = 1:options.maxitr
        V = potential(position);
        g = g_potential(position);

%         fprintf('Itr %04d: %.6e\n', itr, V);

        step = 0.1;
        while step > 1e-6 && potential(position - step*g) > V - options.beta*sum(sum(g.^2))*step
            step = step * options.alpha;
        end
        position = position - step*g;
%         eG(itr) = norm(reshape(g, 2*clusterN, 1));
    end
end