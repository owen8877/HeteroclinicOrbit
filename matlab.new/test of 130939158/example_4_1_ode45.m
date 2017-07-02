B = [[-1 4]; [-4 -1]];
% B = [[-1 0.5]; [-0.5 -1]];
p1 = [4; 1];
x0 = [0; 0];
C = sqrt(2);
h = 1e-3;
sSpan = 0:0.05:1;
xEnd = [];

drawInterval = 5;

%%
xSolution = ode45(@(s, x) [1/sqrt(2); -1/sqrt(2)], [1 0], [0; -1]);
% xSolution = ode45(@(s, x) [0; 0], [0 1], [1; 0]);
for round = 1:50
    pSolution = ode45(@(s, p) dpfunc(s, p, xSolution, B, C), [0 1], p1);    
    xSolution = ode45(@(s, x) dxfunc(s, x, pSolution, B, C), [1 0], [0; -1/(round)]);
    
    p_ = deval(pSolution, sSpan);
    x_ = deval(xSolution, sSpan);
    
    if round - fix(round / drawInterval) * drawInterval == 0
        figure
        s1 = subplot(1,2,1);hold(s1, 'on');
        plot(p_(1, end), p_(2, end), 'r.', 'markersize', 50);
        plot(p_(1, :), p_(2, :));

        s2 = subplot(1,2,2);hold(s2, 'on');
        plot(x_(1, end), x_(2, end), 'r.', 'markersize', 50);
        plot(x_(1, :), x_(2, :));
    end
    
    C = 0;
    for i = 1:(size(sSpan, 2)-1)
        C = C + norm(x_(:, i) - x_(:, i+1));
    end
    fprintf('Round %d, C is %f\n', round, C);
%     xEnd(:, round) = deval(xSolution, 1);
end

function dp = dpfunc(s, p, xSolution, B, C)
    x = deval(xSolution, s);
    b = - B * x;
    dp = C / norm(b) * B' * p;
end

function dx = dxfunc(s, x, pSolution, B, C)
    p = deval(pSolution, s);
    b = - B * x;
    dx = C / norm(b) * (p + b);
end