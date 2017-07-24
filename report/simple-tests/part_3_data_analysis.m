clear
clc

depthRange = 3:4;
B = [[-1 4]; [-4 -1]];
load('high-order-equation-B.mat');

%% Error in s/t param
pError = [];
for depth = depthRange
    load(['s-iteration-' num2str(depth) '.mat']);
    pError(depth) = norm(p1Collection(:, end));
end
figure
hold on;
plot(depthRange, log10(pError(depthRange)));
plot(depthRange, log10(pError(depthRange)), 'ro');

%% Error in phase space; notice that $x$ and $p$ have simple relations
phError = [];
for depth = depthRange
    load(['s-iteration-' num2str(depth) '.mat']);
    
    phErrorArray = zeros(1+10^depth, 1);
    for i = 1:(1+10^depth)
        phErrorArray(i) = norm(pSolution(:, i) - (B + B') * xSolution(:, i));
    end
    phError(depth) = max(phErrorArray);
end
figure
hold on;
plot(depthRange, log10(phError(depthRange)));
plot(depthRange, log10(phError(depthRange)), 'ro');

%% Iterations used
iterations = [];
for depth = depthRange
    load(['s-iteration-' num2str(depth) '.mat']);
    
    xinL = 1;
    while true
        count = size(find(xInitInfo == 2^(-xinL+1)), 2);
        if count == 0
            break
        end
        iterations(depth, xinL) = count;
        xinL = xinL + 1;
    end
end
% figure
% hold on
% for depth = depthRange
%     plot(1:size(iterations, 2), iterations(depth, :), 'DisplayName', num2str(depth));
% end
% hold off
% legend show;

fprintf('\\begin{table}[]\n');
fprintf('\\centering\n');
fprintf('\\caption{Iterations used}\n');
fprintf('\\label{tab:p3-iterations}\n');
fprintf('\\begin{tabular}{c|cccc}\n');
fprintf('\\hline\n');
fprintf('\\multirow{2}{*}{log2($x^{(k)}(1)$)} & \\multicolumn{%d}{c}{Resolution}\\\\\n', size(depthRange, 2));
fprintf('& $10^%d$', depthRange);
fprintf('\\\\\n');
fprintf('\\hline\n');
for i = 1:size(iterations, 2)
    fprintf('$%d$', 1-i);
    fprintf('& %d', iterations(depthRange, i));
    fprintf('\\\\\n');
end
fprintf('\\hline\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{table}\n');