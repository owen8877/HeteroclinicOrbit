% clear
clc
param = getDefaultParam();

%%
resolution = 1e3;
iteration = 0;
iterationLimit = 1;

lStable = [-1; -1];
rStable = [0; 0];

lHess = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4), param), [lStable; 0; 0]);
rHess = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4), param), [rStable; 0; 0]);

%%
[V, D] = eig(lHess);
[Vt, Dt] = eig(rHess);

%% search
scale = 1 / resolution;
ldelta = V(:, 4);
rdelta = Vt(:, 4);

solution = simpleSymplecticSearch(...
    @(~, v, ~) dH(v(1:2), v(3:4), param), ...
    [lStable; 0; 0] + ldelta * scale, [rStable; 0; 0], 1/resolution, 0, -Vt(:, 4));
figure
plotHelper(subplot(2, 2, 1), solution([1, 3], :), 'x', 'px', [[lStable(1) rStable(1)]; [0 0]]);
plotHelper(subplot(2, 2, 2), solution([2, 4], :), 'y', 'py', [[lStable(2) rStable(2)]; [0 0]]);
plotHelper(subplot(2, 2, 3), solution([1, 2], :), 'x', 'y', [lStable rStable]);
plotHelper(subplot(2, 2, 4), solution([3, 4], :), 'px', 'py', [0; 0]);

H = zeros(1, resolution+1);
for i = 1:resolution+1
    H(i) = Hfunc(solution(1:2, i), solution(3:4, i), param);
end

figure
plot(abs(H))

%% itr 1
gridN = size(solution, 2) - 1;
arcSpan = 0;
for i = 1:gridN
    arcSpan = arcSpan + norm(solution(:, i) - solution(:, i+1));
end

solution2 = simpleSymplecticODESolver(...
    @(~, v, u) dH(v(1:2), u(3:4), param), [0, -arcSpan], ...
    [rStable; 0; 0] - rdelta * scale, solution, gridN);
figure
plotHelper(subplot(2, 2, 1), solution2([1, 3], :), 'x', 'px', [[lStable(1) rStable(1)]; [0 0]]);
plotHelper(subplot(2, 2, 2), solution2([2, 4], :), 'y', 'py', [[lStable(2) rStable(2)]; [0 0]]);
plotHelper(subplot(2, 2, 3), solution2([1, 2], :), 'x', 'y', [lStable rStable]);
plotHelper(subplot(2, 2, 4), solution2([3, 4], :), 'px', 'py', [0; 0]);

H = zeros(1, resolution+1);
for i = 1:resolution+1
    H(i) = Hfunc(solution2(1:2, i), solution2(3:4, i), param);
end

figure
plot(abs(H))

%% itr 2
arcSpan = 0;
for i = 1:gridN
    arcSpan = arcSpan + norm(solution(:, i) - solution(:, i+1));
end

solution3 = simpleSymplecticODESolver(...
    @(~, v, u) dH(u(1:2), v(3:4), param), [0, arcSpan], ...
    [lStable; 0; 0] + ldelta * scale, solution, gridN);
figure
plotHelper(subplot(2, 2, 1), solution3([1, 3], :), 'x', 'px', [[lStable(1) rStable(1)]; [0 0]]);
plotHelper(subplot(2, 2, 2), solution3([2, 4], :), 'y', 'py', [[lStable(2) rStable(2)]; [0 0]]);
plotHelper(subplot(2, 2, 3), solution3([1, 2], :), 'x', 'y', [lStable rStable]);
plotHelper(subplot(2, 2, 4), solution3([3, 4], :), 'px', 'py', [0; 0]);

H = zeros(1, resolution+1);
for i = 1:resolution+1
    H(i) = Hfunc(solution3(1:2, i), solution3(3:4, i), param);
end

figure
plot(abs(H))

%%
% function dxdsv = dxds(x, p, param, C)
%     dpHv = dpHn(x, p);
%     H = Hfunc(x, p);
%     dxdsv = dpHv / sqrt(norm(dpHv, 2)^2 - param.lambda * H) * C;
% end
% 
% function dpdsv = dpds(x, p, param, C)
%     dpHv = dpHn(x, p);
%     dxHv = dxHn(x, p);
%     H = Hfunc(x, p);
%     dpdsv = - dxHv / sqrt(norm(dpHv, 2)^2 - param.lambda * H) * C;
% end

function dHv = dH(x, p, param)
    dxH = dxHnum(x, p, param);
    dpH = dpHnum(x, p, param);
    dHv = [dpH; -dxH];
    dHv = dHv / norm(dHv);
end

function dHv = dHsimple(x, p, param)
    dxH = dxHnum(x, p, param);
    dpH = dpHnum(x, p, param);
    dHv = [dpH; -dxH];
end

function dxHv = dxHnum(x, p, param)
    h = 1e-10;
    dxHv = [ ...
        (Hfunc(x+[h; 0], p, param) - Hfunc(x-[h; 0], p, param)) / (2*h); ...
        (Hfunc(x+[0; h], p, param) - Hfunc(x-[0; h], p, param)) / (2*h) ...
        ];
end

function dpHv = dpHnum(x, p, param)
    h = 1e-10;
    dpHv = [ ...
        (Hfunc(x, p+[h; 0], param) - Hfunc(x, p-[h; 0], param)) / (2*h); ...
        (Hfunc(x, p+[0; h], param) - Hfunc(x, p-[0; h], param)) / (2*h) ...
        ];
end

function H = Hfunc(x, p, param)
    H = (p(1)^2 + p(2)^2) / 2 + p(1)*(x(1)-x(1).^3) + p(2)*(x(2)-x(2).^3) + ((x(2)^2-1)*(x(1)^2) - param.L * (x(1)^2-1)*(x(2)^2))^2;
end

%%
function param = getDefaultParam()
    param = struct();
    param.L = 7;
end