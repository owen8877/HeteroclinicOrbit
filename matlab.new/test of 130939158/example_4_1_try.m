clear;
clc;

B = [[-1 4]; [-4 -1]];
p0 = [4; 1];
x1 = [0; 0];

xSolution = struct();
xSolution.t = [0 1];
xSolution.y = [[1; 0] [0; 1]];

C = sqrt(2);

drawInterval = nan;
rounds = nan;

p1Collection = [];

xTErr = [];
pTErr = [];

iteration = 0;
resolution = 1e5;
criterion = 1e0 / resolution;
previousP1 = [];
%% Iteration
while true
    pSolution = myODESolver(@(s, p, x) pDfunc(s, p, x, B, C), [0 1], p0, xSolution, resolution);
    delta = 1e-2 / iteration;
    xSolution = myODESolver(@(s, x, p) xDfunc(s, x, p, B, C), [1 0], x1 + [delta; 0], pSolution, resolution);

    iteration = iteration + 1;
    if iteration - fix(iteration / drawInterval) * drawInterval == 0
        figure
        title(['Iteration ' num2str(iteration)])
        subplot(1,2,1);
        plot(pSolution.y(1, :), pSolution.y(2, :));
        subplot(1,2,2);
        plot(xSolution.y(1, :), xSolution.y(2, :));
    end
    
    %% verify the solution
%     ds = 1e-3;
%     s = 0:ds:1;
% %     x_ = interp1(xSolution.t', xSolution.y', s, 'linear', 'extrap')';
% %     p_ = interp1(pSolution.t', pSolution.y', s, 'linear', 'extrap')';
%     x_ = zeros(size(s, 2), size(xSolution.y, 1));
%     hint = 1;
%     for i = 1:size(s, 2)
%         [x_(i, :), hint] = quickInterpolate(xSolution.t', xSolution.y', s(i), hint);
%     end
%     x_ = x_';
%     
%     p_ = zeros(size(s, 2), size(pSolution.y, 1));
%     hint = 1;
%     for i = 1:size(s, 2)
%         [p_(i, :), hint] = quickInterpolate(pSolution.t', pSolution.y', s(i), hint);
%     end
%     p_ = p_';
%     
%     tArray = zeros(1, size(s, 2));
%     for index = 1:size(s, 2)-1
%         x__ = x_(:, index);
%         p__ = p_(:, index);
%         b = -B * x__;
%         lambda = norm(p__ + b) / C;
%         tArray(index+1) = tArray(index) + 1/lambda * ds;
%     end
%     
%     Rp = zeros(2, size(s, 2));
%     Rx = zeros(2, size(s, 2));
%     pErrorArray = zeros(1, size(s, 2));
%     xErrorArray = zeros(1, size(s, 2));
%     for index = 1:size(s, 2)
%         Rp(:, index) = expm(B' * tArray(index)) * p1;
%         Rx(:, index) = (B + B') \ Rp(:, index);
%         pErrorArray(index) = norm(Rp(:, index) - p_(:, index));
%         xErrorArray(index) = norm(Rx(:, index) - x_(:, index));
%     end
%     
%     pTErr(round) = max(pErrorArray);
%     xTErr(round) = max(xErrorArray);
    
    %% update C
    C = 0;
    for i = 1:(size(xSolution.y, 2) - 1)
        C = C + norm(xSolution.y(:, i) - xSolution.y(:, i+1));
    end
    clc
    fprintf('Iteration %d completed.\n', iteration);
    
    if iteration > rounds
        break
    end
    
    %% watch if p(1) converges
    if iteration ~= 1
        p1 = pSolution.y(:, end);
        endDistance = norm(p1 - previousP1);
        if endDistance < criterion
            break
        end
        fprintf('probable error in p1: %e', endDistance);
    end
    previousP1 = pSolution.y(:, end);
    p1Collection(:, iteration) = previousP1;
end

% figure
% plot(1:rounds, log10(xTErr));
% title('xTErr');
% figure
% plot(1:rounds, log10(pTErr));
% title('pTErr');

function pD = pDfunc(~, p, x, B, C)
    b = - B * x;
    bM = norm(b);
    pD = C / bM * B' * p;
end

function xD = xDfunc(~, x, p, B, C)
    b = - B * x;
    bM = norm(b);
    xD = C / bM * (p + b);
end