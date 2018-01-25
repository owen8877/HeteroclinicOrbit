clear
% clc
param = getDefaultParam();

%% we first find path from p1 to p2 above x-axis
xLeft = -1;
xRight = 0;
resolution = 1e3;
xSolution = linspace(xLeft, xRight, resolution+1);
C = xRight - xLeft;
iteration = 0;
iterationLimit = 200;
drawInterval = 200;

pErr = zeros(iterationLimit, 1);
HErr = zeros(iterationLimit, 1);
while iteration < iterationLimit
    pSolution = naiveODESolver(@(s, p, x) dpds(x, p, param, C), [1 0], 0, xSolution, resolution, ...
        @(p, x) characterDirection(x, p));
%     xSolution = naiveODESolver(@(s, x, p) dxds(x, p, param, C), [0 1], xLeft, pSolution, resolution, ...
%         @(x, p) C);
    iteration = iteration + 1;
    
    C = 0;
    for i = 1:resolution
        C = C + norm(xSolution(:, i) - xSolution(:, i+1));
    end
    
    pPrecise = preciseSolution(xSolution);
    pErr(iteration) = max(abs(pSolution - pPrecise));
    H = zeros(resolution+1, 1);
    for i = 1:numel(H)
        H(i) = Hfunc(xSolution(:, i), pSolution(:, i));
    end
    HErr(iteration) = max(abs(H));
    
    if drawInterval > 0 && mod(iteration, drawInterval) == 0
        figure
        hold on
        plot(xSolution, pSolution);
        plot(xSolution, xSolution .* (xSolution .^ 2 - 1) * 2, 'r');
        
%         figure
%         plot(H);
    end
    
    if true
        clc
        fprintf('Iteration %d / %d \n', iteration, iterationLimit);
    end
end

figure
plot(pErr)
figure
plot(HErr)

figure
hold on
plot(xSolution, pSolution);
plot(xSolution, xSolution .* (xSolution .^ 2 - 1) * 2, 'r')

%%

function pPrecise = preciseSolution(xSpan)
    pPrecise = xSpan * 0;
    for i = 1:size(xSpan, 2)
        pPrecise(i) = 2 * xSpan(i)^3 - 2 * xSpan(i);
    end
end

function slope = characterDirection(x, p)
    slope = -2 * dpdxHn(x, p) / dpdpHn(x, p);
end

function dpdxHv = dpdxHn(x, p)
    h = 1e-8;
    dpdxHv = (dxH(x, p+h) - dxH(x, p-h)) / 2*h;
end

function dpdpHv = dpdpHn(x, p)
    h = 1e-8;
    dpdpHv = (dpH(x, p+h) - dpH(x, p-h)) / 2*h;
end

function dxdsv = dxds(x, p, param, C)
    dpHrv = dpH(x, p);
    H = Hfunc(x, p);
    dxdsv = dpHrv / sqrt(dpHrv^2) * C;
end

function dpdsv = dpds(x, p, param, C)
    dpHv = dpH(x, p);
    dxHv = dxH(x, p);
    H = Hfunc(x, p);
    dpdsv = - dxHv / sqrt(dpHv^2 - param.lambda * H) * C;
end

function dpHv = dpH(x, p)
    dpHv = p + x - x^3;
end

function dxHv = dxH(x, p)
    dxHv = p * (1 - 3 * x^2);
end

function Hv = Hfunc(x, p)
    Hv = p * p / 2 + p * (x - x^3);
end

function param = getDefaultParam()
    param = struct();
    param.lambda = 2;
end