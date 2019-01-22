clear; %clc
addpath util

clusterN = 7;

global vPad hPad
vPad = repmat(1:clusterN, clusterN, 1);
hPad = vPad';

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
fixA = 0;
positionB = [ ...
    0, 0; ...
    s62*1.5, -s62*s3d2; ...
    s62, 0; ...
    s62/2, s62*s3d2; ...
    -s62, 0; ...
    -s62/2, -s62*s3d2; ...
    s62/2, -s62*s3d2; ...
]'; % State B
fixB = 0;
positionC = [ ...
    0, 0; ...
    s62*2, 0; ...
    s62, 0; ...
    s62/2, s62*s3d2; ...
    -s62/2, -s62*s3d2; ...
    s62/2, -s62*s3d2; ...
    s62*1.5, -s62*s3d2; ...
]'; % State C
positionD = [ ...
    -s62/2, s62*s3d2; ...
    s62, 0; ...
    0, 0; ...
    s62/2, s62*s3d2; ...
    -s62, 0; ...
    -s62/2, -s62*s3d2; ...
    s62/2, -s62*s3d2; ...
]'; % State D
positionA = reshape(centerize(positionA), 2*clusterN, 1);
positionB = reshape(centerize(positionB), 2*clusterN, 1);
positionC = reshape(centerize(positionC), 2*clusterN, 1);
positionD = reshape(centerize(positionD), 2*clusterN, 1);

options.maxitr = 100;
options.beta = 0.5;
options.alpha = 0.5;
positionA = initialSolver(positionA, options);
positionB = initialSolver(positionB, options);
positionC = initialSolver(positionC, options);
positionD = initialSolver(positionD, options);

%% React-Diff Equation
h = 1/16;
grid = (-12+h/2):h:(12-h/2);
points = [positionA positionD]';

if exist('backup.old/matlabAD.mat', 'file')
    load('backup.old/matlabAD.mat');
else
%     sol = initialGenerator(grid, points);
    sol = mergeInit();
end

sol = rotate(sol);

time = 0;
timeStep = 2e-5;
itr = 1;
maxItr = 2000;

for i = 1:4
    figure(i); clf
end

results = zeros(maxItr, 2*clusterN, numel(grid));
tic
errorV = zeros(maxItr/100, 1);
rhsE = zeros(maxItr/100, 1);
while itr <= maxItr
    lsol = laplacian(sol, h);
    updateQ = ReactionTerm(sol) + lsol;
    sol = sol + updateQ * timeStep;
    
    results(itr, :, :) = sol;
    
    time = time + timeStep;
    itr = itr + 1;
    
    if mod(itr, floor(maxItr/5)) == 0
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
        drawnow

        p = derivative(sol, h) - hV;
        hamiltonV = hamiltonianFunc(sol, p);
        figure(3); hold on
        plot(hamiltonV)
        drawnow

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
        drawnow
    end

    if mod(itr, 100) == 0
        p = derivative(sol, h) - hFunc(sol);
        hamiltonV = hamiltonianFunc(sol, p);
        fprintf('Iteration % 6d, error in H %.6f.\n', itr, max(abs(hamiltonV)));
        errorV(itr/100) = max(abs(hamiltonV));
        rhsE(itr/100) = max(max(abs(updateQ)));
    end
end
toc

% save('backup.old/matlabAD.mat', 'sol')

figure(2); legend('1', '2', '3', '4', '5')
figure(3); legend('1', '2', '3', '4', '5')
figure(4); legend('1', '2', '3', '4', '5')

figure(1)
contour(permute(results(:, 1, :), [1 3 2]))

function f = ReactionTerm(q)
%     f = nlFunc(q, fFunc) - nhhFunc(q, fFunc);
    [m, n] = size(q);
    f = zeros(m, n);
    parfor i = 1:n
        f(:, i) = - nhFuncSimp(q(:, i))' * hFuncSimp(q(:, i));
    end
end

function hv = hFunc(q)
    [m, n] = size(q);
    hv = zeros(m, n);
    parfor i = 1:n
        hv(:, i) = g_potential(q(:, i));
    end
end

function hv = hFuncSimp(q)
    hv = g_potential(q);
end

function nhm = nhFunc(q)
    [m, n] = size(q);
    nhm = zeros(m, m, n);
    parfor i = 1:n
        nhm(:, :, i) = h_potential(q(:, i));
    end
end

function nhm = nhFuncSimp(q)
    nhm = h_potential(q);
end

function h = hamiltonianFunc(q, p)
    h = sum(p.^2/2 + p.*hFunc(q));
%     [~, n] = size(q);
%     h = zeros(1, n);
%     for i = 1:n
%         qv = q(:, i);
%         pv = p(:, i);
%         h(:, i) = pv'*pv/2 + pv'*hFunc(qv);
%     end
end

%% Initial Boundary
function position = centerize(position)
    center = mean(position, 2);
    position = position - center;
end

function V = potential(position)
    n = numel(position)/2;
    d = distance(position);
    
    halfd = reshape(triu(d, 1), n^2, 1);
    halfd = halfd(halfd ~= 0);
    
    V = sum(halfd.^-12-halfd.^-6);
end

function g = g_potential(position)
    global vPad hPad
    x = position(1:2:end); y = position(2:2:end);
    n = numel(position)/2;
    d = distance(position);
    
    g = zeros(2*n, 1);
    fddr = 6./d.^8 - 12./d.^14;
    fddr(1:1+n:end) = 0;
    dxM = x(vPad)-x(hPad); dyM = y(vPad)-y(hPad);
    for i = 1:n
        g(2*i-1) = fddr(i, :) * dxM(i, :)';
        g(2*i) = fddr(i, :) * dyM(i, :)';
    end
end

function h = h_potential(position)
    global vPad hPad
    x = position(1:2:end); y = position(2:2:end);
    n = numel(position)/2;
    d = distance(position);
    
    fddr = 6./d.^8 - 12./d.^14;
    fddr(1:1+n:end) = 0;
    dxM = x(vPad)-x(hPad); dyM = y(vPad)-y(hPad);
    
    fddrddr = 48./d.^10 - 168./d.^16;
    fddrddr(1:1+n:end) = 0;
    hxx = -(fddrddr .* dxM.^2 + fddr);
    hxx(1:1+n:end) = 0;
    hxx(1:1+n:end) = -sum(hxx, 2);
    hxy = -(fddrddr .* dxM .* dyM);
    hxy(1:1+n:end) = 0;
    hxy(1:1+n:end) = -sum(hxy, 2);
    hyx = -(fddrddr .* dxM .* dyM);
    hyx(1:1+n:end) = 0;
    hyx(1:1+n:end) = -sum(hyx, 2);
    hyy = -(fddrddr .* dyM.^2 + fddr);
    hyy(1:1+n:end) = 0;
    hyy(1:1+n:end) = -sum(hyy, 2);
    
    h = zeros(2*n);
    h(1:2:end, 1:2:end) = hxx;
    h(1:2:end, 2:2:end) = hxy;
    h(2:2:end, 1:2:end) = hyx;
    h(2:2:end, 2:2:end) = hyy;
end

function d = distance(position) % x1 y1 x2 y2 ...
    global vPad hPad
    x = position(1:2:end); y = position(2:2:end);
    d = sqrt((x(vPad)-x(hPad)).^2 + (y(vPad)-y(hPad)).^2);
end

function position = initialSolver(position, options)
    for itr = 1:options.maxitr
        V = potential(position);
        g = g_potential(position);

%         fprintf('Itr %04d: %.6e\n', itr, V);

        step = 0.1;
        while step > 1e-6 && potential(position - step*g) > V - options.beta*(g'*g)*step
            step = step * options.alpha;
        end
        position = position - step*g;
%         eG(itr) = norm(reshape(g, 2*clusterN, 1));
    end
end