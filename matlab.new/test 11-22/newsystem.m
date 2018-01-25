clc; clear;

% init
q1Stable = -1/sqrt(2);
q1Saddle = 0;
lVal = [q1Stable; 0; 0; 0];
rVal = [q1Saddle; 0; 0; 0];
rVal2 = rVal;
resolution = 1e4;

lHessian = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4)), lVal);
[lV, lD] = eig(lHessian);
rHessian = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4)), rVal);
[rV, rD] = eig(rHessian);

%% Search part
delta = real(lV(:, 2)) / resolution;

V = lV; V(:, [1 2]) = lV(:, [2 1]);
solution = simpleSymplecticSearch(@(~, v, ~) dHn(v(1:2), v(3:4)), ...
    lVal+delta, rVal, 1/resolution, 0, rV(:, 1));

figure
plotHelper(subplot(2, 2, 1), solution([1 3], :), 'q1', 'p1', [lVal([1 3]) rVal([1 3]) rVal2([1 3], :)]);
plotHelper(subplot(2, 2, 2), solution([2 4], :), 'q2', 'p2', [lVal([2 4]) rVal([2 4]) rVal2([2 4], :)]);
plotHelper(subplot(2, 2, 3), solution([1 2], :), 'q1', 'q2', [lVal([1 2]) rVal([1 2]) rVal2([1 2], :)]);
plotHelper(subplot(2, 2, 4), solution([3 4], :), 'p1', 'p2', [lVal([3 4]) rVal([3 4]) rVal2([3 4], :)]);

H = zeros(1, size(solution, 2));
for i = 1:size(solution, 2)
    H(i) = Hfunc(solution(1:2, i), solution(3:4, i));
end

figure
plot(abs(H))

%% Search from right
% alpha = 0.5;
% delta = (alpha * real(rV(:, 3)) + (1-alpha) * real(rV(:, 4))) / resolution;
% V = lV; V(:, [1 4]) = lV(:, [1 4]);
% solution = simpleSymplecticODESolver(@(~, v, ~) dH(v(1:2), v(3:4)), ...
%     [0, 5], rVal+delta, [], resolution);
% % solution = simpleSymplecticSearch(@(~, v, ~) dH(v(1:2), v(3:4)), ...
% %     rVal+delta, lVal, 1/resolution, 1, V);
% 
% figure
% plotHelper(subplot(2, 2, 1), solution([1 3], :), 'q1', 'p1', [lVal([1 3]) rVal([1 3]) rVal2([1 3], :)]);
% plotHelper(subplot(2, 2, 2), solution([2 4], :), 'q2', 'p2', [lVal([2 4]) rVal([2 4]) rVal2([2 4], :)]);
% plotHelper(subplot(2, 2, 3), solution([1 2], :), 'q1', 'q2', [lVal([1 2]) rVal([1 2]) rVal2([1 2], :)]);
% plotHelper(subplot(2, 2, 4), solution([3 4], :), 'p1', 'p2', [lVal([3 4]) rVal([3 4]) rVal2([3 4], :)]);
% 
% H = zeros(1, size(solution, 2));
% for i = 1:size(solution, 2)
%     H(i) = Hfunc(solution(1:2, i), solution(3:4, i));
% end
% 
% figure
% plot(abs(H))

%% Reverse search
% solution = simpleSymplecticODESolver(@(~, v, ~) dH(v(1:2), v(3:4)), ...
%     [0, 0.6], lVal + delta, [], resolution);
% alpha = 1;
% delta = (alpha * real(rV(:, 1)) + (1-alpha) * real(rV(:, 2))) / resolution;
% V = lV; V(:, [1 4]) = lV(:, [1 4]);
% solution = simpleSymplecticSearch(@(~, v, ~) dHn(v(1:2), v(3:4)), ...
%     rVal+delta, lVal, - 1/resolution, 1, V);
% 
% figure
% plotHelper(subplot(2, 2, 1), solution([1 3], :), 'q1', 'p1', [lVal([1 3]) rVal([1 3]) rVal2([1 3], :)]);
% plotHelper(subplot(2, 2, 2), solution([2 4], :), 'q2', 'p2', [lVal([2 4]) rVal([2 4]) rVal2([2 4], :)]);
% plotHelper(subplot(2, 2, 3), solution([1 2], :), 'q1', 'q2', [lVal([1 2]) rVal([1 2]) rVal2([1 2], :)]);
% plotHelper(subplot(2, 2, 4), solution([3 4], :), 'p1', 'p2', [lVal([3 4]) rVal([3 4]) rVal2([3 4], :)]);
% 
% H = zeros(1, size(solution, 2));
% for i = 1:size(solution, 2)
%     H(i) = Hfunc(solution(1:2, i), solution(3:4, i));
% end
% 
% figure
% plot(abs(H))

%% 1st itr: know q, solve p
% itrResolution = size(solution, 2) - 1;
% bSolution = simpleSymplecticODESolver(@(~, x, z) Hfunc(z, x), [1 0], ...
%     rVal(3:4), solution(1:2, :), itrResolution);
% 
% figure
% plotHelper(subplot(2, 2, 1), [solution(1, :); bSolution(1, :)], 'q1', 'p1', [lVal([1 3]) rVal([1 3]) rVal2([1 3], :)]);
% plotHelper(subplot(2, 2, 2), [solution(2, :); bSolution(2, :)], 'q2', 'p2', [lVal([2 4]) rVal([2 4]) rVal2([2 4], :)]);
% plotHelper(subplot(2, 2, 3), [solution(1, :); solution(2, :)], 'q1', 'q2', [lVal([1 2]) rVal([1 2]) rVal2([1 2], :)]);
% plotHelper(subplot(2, 2, 4), [bSolution(1, :); bSolution(2, :)], 'p1', 'p2', [lVal([3 4]) rVal([3 4]) rVal2([3 4], :)]);
% 
% for i = 1:size(solution, 2)
%     H(i) = Hfunc(solution(1:2, i), bSolution(:, i));
% end
% 
% figure
% plot(abs(H))

%%
function dHv = dH(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
end

function dHnv = dHn(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
    dHnv = dHv / norm(dHv);
end

function Hpv = Hpfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    Hpv = [q2; p2/2+df(q1)*q2];
end

function Hqv = Hqfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    Hqv = [p2*q2*ddf(q1)+2*(q2-f(q1))*df(q1); p1+p2*df(q1)-2*(q2-f(q1))];
end

function H = Hfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    H = p1*q2 + p2*(p2/4+df(q1)*q2) - (q2-f(q1))^2;
end

function fv = f(x)
    % fv = -(5*x^4 + 8*x^3 + 3*x^2);
    % fv = 12*x^3 + 12*x^2;
    % fv = 2*x^3 - x/2;
    fv = 4*x^3 - 2*x;
end

function dfv = df(x)
    % dfv = -(20*x^3 + 24*x^2 + 6*x);
    % dfv = 36*x^2 + 24;
    % dfv = 6*x^2 - 1/2;
    dfv = 12*x^2 - 2;
end

function ddfv = ddf(x)
    % ddfv = -(60*x^2 + 48*x + 6);
    % ddfv = 72*x + 24;
    % ddfv = 12*x;
    ddfv = 24*x;
end