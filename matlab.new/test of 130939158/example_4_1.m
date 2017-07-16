clear;
clc;

B = [[-1 4]; [-4 -1]];
% B = [[-1 0.5]; [-0.5 -1]];
p1 = [4; 1];
x0 = [0; 0];

xSolution = struct();
xSolution.t = [0 1];
xSolution.y = [[1; 0] [0; 1]];

C = sqrt(2);
% h = 1e-3;

xEnd = [];
drawInterval = 20;
rounds = 100;

xTErr = zeros(1, rounds);
pTErr = zeros(1, rounds);

%% Iteration
for round = 1:rounds
    pSolution = myODESolver(@(s, p, x) pDfunc(s, p, x, B, C), [0 1], p1, xSolution, round);
    xSolution = myODESolver(@(s, x, p) xDfunc(s, x, p, B, C), [1 0], [0; 1/round], pSolution, round);

    if round - fix(round / drawInterval) * drawInterval == 0
        figure
        title(['Round ' num2str(round)])
        subplot(1,2,1);
        plot(pSolution.y(1, :), pSolution.y(2, :));
        subplot(1,2,2);
        plot(xSolution.y(1, :), xSolution.y(2, :));
    end
    
    %% verify the solution
    ds = 1e-3;
    s = 0:ds:1;
%     x_ = interp1(xSolution.t', xSolution.y', s, 'linear', 'extrap')';
%     p_ = interp1(pSolution.t', pSolution.y', s, 'linear', 'extrap')';
    x_ = zeros(size(s, 2), size(xSolution.y, 1));
    hint = 1;
    for i = 1:size(s, 2)
        [x_(i, :), hint] = quickInterpolate(xSolution.t', xSolution.y', s(i), hint);
    end
    x_ = x_';
    
    p_ = zeros(size(s, 2), size(pSolution.y, 1));
    hint = 1;
    for i = 1:size(s, 2)
        [p_(i, :), hint] = quickInterpolate(pSolution.t', pSolution.y', s(i), hint);
    end
    p_ = p_';
    
    tArray = zeros(1, size(s, 2));
    for index = 1:size(s, 2)-1
        x__ = x_(:, index);
        p__ = p_(:, index);
        b = -B * x__;
        lambda = norm(p__ + b) / C;
        tArray(index+1) = tArray(index) + 1/lambda * ds;
    end
    
    Rp = zeros(2, size(s, 2));
    Rx = zeros(2, size(s, 2));
    pErrorArray = zeros(1, size(s, 2));
    xErrorArray = zeros(1, size(s, 2));
    for index = 1:size(s, 2)
        Rp(:, index) = expm(B' * tArray(index)) * p1;
        Rx(:, index) = (B + B') \ Rp(:, index);
        pErrorArray(index) = norm(Rp(:, index) - p_(:, index));
        xErrorArray(index) = norm(Rx(:, index) - x_(:, index));
    end
    
    pTErr(round) = max(pErrorArray);
    xTErr(round) = max(xErrorArray);
    
    %%
    C = 0;
    for i = 1:(size(xSolution.y, 2) - 1)
        C = C + norm(xSolution.y(:, i) - xSolution.y(:, i+1));
    end
    fprintf('Round %d completed.\n', round);
    
    xEnd(:, round) = xSolution.y(:, end);
end

figure
plot(1:rounds, log10(xTErr));
title('xTErr');
figure
plot(1:rounds, log10(pTErr));
title('pTErr');

function pD = pDfunc(s, p, x, B, C)
    b = - B * x;
    bM = norm(b);
    pD = C / bM * B' * p;
end

function xD = xDfunc(s, x, p, B, C)
    b = - B * x;
    bM = norm(b);
    xD = C / bM * (p + b);
end