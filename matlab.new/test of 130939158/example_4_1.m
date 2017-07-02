B = [[-1 4]; [-4 -1]];
% B = [[-1 0.5]; [-0.5 -1]];
p1 = [4; 1];
x0 = [0; 0];
sx = [0 1];
xs = [[0; -1] [-1; 0]];
C = sqrt(2);
h = 1e-3;

xEnd = [];
drawInterval = 10;

%% Iteration
for round = 1:100
    sp = [];
    ps = [p1];

    index = 1;
    s = 0;
    while true
        x = interp1(sx', xs', s, 'linear', 'extrap')';
        b = - B * x;
        bM = norm(b);
        dp = C / bM * B' * ps(:, index);
        ps(:, index + 1) = ps(:, index) + dp * h;
        sp(index) = s;

        s = s + h;
        index = index + 1;
        if s >= 1
            break;
        end
    end
    sp(index) = s;

    %%
    sx = [];
    xs = [0; -1/round];

    index = 1;
    s = 1;
    while true
        p = interp1(sp', ps', s, 'linear', 'extrap')';
        b = - B * xs(:, index);
        bM = norm(b);
        dx = C / bM * (p + b);
        xs(:, index + 1) = xs(:, index) - dx * h;
        sx(index) = s;

        s = s - h;
        index = index + 1;
        if s <= 0
            break;
        end
    end
    sx(index) = s;

    if round - fix(round / drawInterval) * drawInterval == 0
        figure
        title(['Round ' num2str(round)])
        subplot(1,2,1);
        plot(ps(1, :), ps(2, :));
        subplot(1,2,2);
        plot(xs(1, :), xs(2, :));
    end
    
    %%
    C = 0;
    for i = 1:(size(sx, 2)-1)
        C = C + norm(xs(:, i) - xs(:, i+1));
    end
    fprintf('Round %d completed.\n', round);
    
    xEnd(:, round) = interp1(sx', xs', 1, 'linear', 'extrap')';
end

%%
% C = 1;
% 
% [~, p] = ode45(@(s, y) dyfunc(s, y, B, C), [1 0], p0);
% 
% function dy = dyfunc(s, y, B, C)
%     x = [s; 0];
%     b = - B * x;
%     dy = C / norm(b) * B * y;
% end
