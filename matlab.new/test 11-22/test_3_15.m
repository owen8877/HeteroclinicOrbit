clc; clear global; clear

%%
% Preparations
addpath testGen
global Hqfunc Hpfunc Hfunc
[lVal, rVal, Hqfunc, Hpfunc, Hfunc] = selfTestGen(false);

resolution = 1000;
leftSide = false;
gridSearch = true;
search618Step = 1;

lHessian = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4)), lVal);
[lV, lD] = eigs(lHessian, size(lHessian, 1), 'lr');
rHessian = notOrdinaryHessian(@(v) Hfunc(v(1:2), v(3:4)), rVal);
[rV, rD] = eigs(rHessian, size(rHessian, 1), 'lr');

disp([lD rD]);
disp([lV rV]);

%%
% The following code gives a set of possible initial angle between eigen 
% directions

fais = linspace(0, 2*pi, 100)';
fais = fais(2:end);
distance = [];
notice = [];
if gridSearch
    for fai = fais'
        [thisDistance, solution] = simplexSearchWrapper(fai, leftSide, lVal, rVal, ...
            lV, rV, lD, rD, @dHn, resolution);
        distance = [distance; thisDistance];
        fprintf('Searching %.6e with distance %.4e\n', fai, thisDistance);
        if thisDistance < 1
            notice = [notice; fai];
            fprintf('Found one %.6e\n', fai);
        end
    end
    figure
    plot(fais, distance);
end

%%
% The following code searches for the minimal value

% initial = notice(1);
% initial = 0.4649557127; % for self-made cases on the right, L=3
initial = 0.3769; % for self-made cases on the right, L=5
% initial = 0.3896; % with 0.2011, 0.3896, 0.5027, 0.8168, 1.2315, 4.9763
% a possible homoclinic trajectory
% initials = [
%     0.0754
%     0.1634
%     0.4021
%     0.5404
%     0.6660
%     0.9173
%     1.0053
%     1.3320
%     1.7216
%     1.8347
%     1.8975
%     2.1865
%     2.4379
%     2.6389
%     2.8400
%     2.9657
%     ];
% initial = 4.1552; % for prob cases
% initial = 0;
% initial = 2; % for self-made cases on the left, L=3
step = pi/1000;

[initialDis, solution] = simplexSearchWrapper(initial, leftSide, lVal, rVal, ...
    lV, rV, lD, rD, @dHn, resolution);
fprintf('Init dis %.4e\n', initialDis);

for leftMulti = 1:20
    fai = initial - (leftMulti+1) * step;
    [dis, ~] = simplexSearchWrapper(fai, leftSide, lVal, rVal, ...
        lV, rV, lD, rD, @dHn, resolution);
    fprintf('LeftMulti %d distance %.4e\n', leftMulti, dis);
    if dis > initialDis
        break
    end
end
for i = 1:10
    left = initial - (leftMulti+2^(-i)) * step;
    [dis, ~] = simplexSearchWrapper(left, leftSide, lVal, rVal, ...
        lV, rV, lD, rD, @dHn, resolution);
    fprintf('Lefti %d distance %.4e\n', i, dis);
    if dis < initialDis
        break
    end
end

for rightMulti = 1:20
    fai = initial + (rightMulti+1) * step;
    [dis, ~] = simplexSearchWrapper(fai, leftSide, lVal, rVal, ...
        lV, rV, lD, rD, @dHn, resolution);
    fprintf('RightMulti %d distance %.4e\n', rightMulti, dis);
    if dis > initialDis
        break
    end
end
for i = 1:10
    right = initial + (rightMulti+2^(-i)) * step;
    [dis, ~] = simplexSearchWrapper(right, leftSide, lVal, rVal, ...
        lV, rV, lD, rD, @dHn, resolution);
    fprintf('Righti %d distance %.4e\n', i, dis);
    if dis < initialDis
        break
    end
end

leftValueKnown = false; rightValueKnown = false;
for step = 1:search618Step
    midleft = right - (right-left) * (sqrt(5)-1)/2;
    midright = left + (right-left) * (sqrt(5)-1)/2;
    if ~leftValueKnown
        [midleftValue, solution] = ...
            simplexSearchWrapper(midleft, leftSide, ...
            lVal, rVal, lV, rV, lD, rD, @dHn, resolution);
    end
    
    if ~rightValueKnown
        [midrightValue, solution] = ...
            simplexSearchWrapper(midright, leftSide, ...
            lVal, rVal, lV, rV, lD, rD, @dHn, resolution);
    end
    
    fprintf('%.6e, %.6e - %.4e, %.4e\n', left, right, midleftValue, midrightValue);
        
    if midleftValue < midrightValue
        right = midright;
        midrightValue = midleftValue;
        rightValueKnown = true;
        leftValueKnown = false;
    else
        left = midleft;
        midleftValue = midrightValue;
        rightValueKnown = false;
        leftValueKnown = true;
    end
end
richPlotHelper(0, solution, lVal, rVal, lVal, Hfunc, @dH, ...
    struct('extraPlot', true));
fprintf('initial %.4e final %.4e\n', initial, left);
save('data/solution.mat', 'solution');

%%
function [distance, solution] = simplexSearchWrapper(fai, leftSide, lVal, rVal, lV, rV, lD, rD, dHn, resolution)
    global Hfunc
    if ~leftSide
        delta = real(rV * (rD \ [cos(fai); sin(fai); 0; 0]));
        delta = delta / norm(delta);
        rAnchor = rVal - delta / resolution;
        % rAnchor = rVal + (rV(:, 1)*sin(fai) + rV(:, 2)*cos(fai)) / resolution;
        lV_ = lV;
        lV_(:, [1 4]) = lV_(:, [4 1]);
        solution = simpleSymplecticSearch(@(~, v, ~) dHn(v(1:2), v(3:4)), ...
            rAnchor, lVal, -1/resolution, 0, real(lV_(:, 1)), struct('output', false));
        solution = fliplr(solution);
        distance = norm(lVal - solution(:, 1));
    else
        delta = real(lV * (lD \ [0; 0; cos(fai); sin(fai)]));
        delta = delta / norm(delta);
        lAnchor = lVal + delta / resolution;
        solution = simpleSymplecticSearch(@(~, v, ~) dHn(v(1:2), v(3:4)), ...
            lAnchor, rVal, 1/resolution, 0, real(rV(:, 1)), struct('output', false));
        distance = norm(rVal - solution(:, end));
    end
%     richPlotHelper(0, solution, lVal, rVal, lVal, Hfunc, @dH, struct());
end

%%
function nv = normalize(v)
    nv = v / norm(v);
end

function dHv = dH(q, p)
    global Hqfunc Hpfunc
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
end

function dHv = nablaH(q, p)
    global Hqfunc Hpfunc
    dHv = [Hqfunc(q, p); Hpfunc(q, p)];
end

function dHnv = dHn(q, p)
    global Hqfunc Hpfunc
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
    dHnv = dHv / norm(dHv);
end

function dHdpnv = dHdpn(q, p)
    global Hqfunc Hpfunc
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
    dHnv = dHv / norm(dHv);
    dHdpnv = dHnv(3:4);
end

function dHdqnv = dHdqn(q, p)
    global Hqfunc Hpfunc
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
    dHnv = dHv / norm(dHv);
    dHdqnv = dHnv(1:2);
end