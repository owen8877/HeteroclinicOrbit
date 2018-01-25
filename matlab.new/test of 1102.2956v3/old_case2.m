clear
% clc
param = getDefaultParam();

%% we first find path from p1 to p2 above x-axis
xLeft = [0; 0];
xRight = [1; 1];
resolution = 1e3;
xSolution = zeros(2, resolution+1);
xSolution(2, :) = linspace(0, 1, resolution+1);
xSolution(1, :) = xSolution(2, :) .^ 1;
C = solutionLength(xSolution);
iteration = 0;
iterationLimit = 5;
drawInterval = 1;

% pErr = zeros(iterationLimit, 1);
% while iteration < iterationLimit
%     iteration = iteration + 1;
%     
%     pSolution = naiveODESolver(@(s, p, x) dpds(x, p, param, C), [1 0], ...
%         [0; 0], xSolution, resolution, @(p, x) charDirInP(x, p, C));
%     xSolution = naiveODESolver(@(s, x, p) dxds(x, p, param, C), [0 1], ...
%         xLeft, pSolution, resolution, @(x, p) charDirInX(x, p, C));
%     
%     C = solutionLength(xSolution);
%     
%     % pErr(iteration) = max(abs(pSolution - pPrecise));
%     
%     if drawInterval > 0 && mod(iteration, drawInterval) == 0
%         figure
%         s1 = subplot(1, 2, 1);
%         hold(s1, 'on');
%         plot(xSolution(1, :), xSolution(2, :));
%         plot(xSolution(1, 1), xSolution(2, 1), 'ro')
%         s2 = subplot(1, 2, 2);
%         hold(s2, 'on');
%         plot(pSolution(1, :), pSolution(2, :));
%         plot(pSolution(1, end), pSolution(2, end), 'ro')
%         
%         H = zeros(1, resolution+1);
%         for i = 1:resolution+1
%             H(i) = Hfunc(xSolution(:, i), pSolution(:, i));
%         end
%         figure
%         plot(H)
%     end
%     
%     if true
%         clc
%         fprintf('Iteration %d / %d \n', iteration, iterationLimit);
%     end
% end

solution = simpleSymplecticODESolver(@(~, v, ~) dH(v(1:2), v(3:4)), [0, 3], ...
    [-0.99; -0.99; 0.04; 0.04], xSolution, resolution);
figure
plot(solution(1, :), solution(3, :));

figure
plot(solution(1, :), solution(2, :));

H = zeros(1, resolution+1);
for i = 1:resolution+1
    H(i) = Hfunc(solution(1:2, i), solution(3:4, i));
end

figure
plot(H)

%%
function slope = characterDirection(x, px)
    slope = -2 * dpdxHrn(x, px) / dpdpHrn(x, px);
end

function charD = charDirInX(x, p, C)
    Hesslike = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4)), [x; p]);
    [V, D] = eig(Hesslike);
    fullCharD = [1; 0; 0; 0];
    charD = fullCharD / norm(fullCharD(1:2)) * C;
    charD = charD(1:2);
end

function charD = charDirInP(x, p, C)
    Hesslike = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4)), [x; p]);
    [V, D] = eig(Hesslike);
    fullCharD = [-1; 1; 6; -6];
    charD = fullCharD / norm(fullCharD(1:2)) * C;
    charD = charD(3:4) / 2;
end

function dxdsv = dxds(x, p, param, C)
    dpHv = dpH(x, p);
    H = Hfunc(x, p);
    dxdsv = dpHv / sqrt(norm(dpHv, 2)^2 - param.lambda * H) * C;
end

function dpdsv = dpds(x, p, param, C)
    dpHv = dpH(x, p);
    dxHv = dxH(x, p);
    H = Hfunc(x, p);
    dpdsv = - dxHv / sqrt(norm(dpHv, 2)^2 - param.lambda * H) * C;
end

function dHv = dH(x, p)
    dxHv = dxH(x, p);
    dpHv = dpH(x, p);
    dHv = [dpHv; -dxHv] / norm([dpHv; -dxHv]);
end

function dpHv = dpH(x, p)
    dpHv = [p(1) + x(1) - x(1)^3; p(2) + x(2) - x(2)^3];
end

function dxHv = dxH(x, p)
    dxHv = [p(1) - 3*x(1)^2*p(1); p(2) - 3*x(2)^2*p(2)];
end

function H = Hfunc(x, p)
    H = (p(1)^2 + p(2)^2) / 2 + p(1)*(x(1)-x(1).^3) + p(2)*(x(2)-x(2).^3);
end

function param = getDefaultParam()
    param = struct();
    param.lambda = 0;
end