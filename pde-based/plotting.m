clear; % clc

% % fig:pde-dir-vs-neu
% dirichlet = load('plot-data/simple-dirichlet.mat');
% neumann = load('plot-data/simple-neumann.mat');
% figure; grid on;
% semilogy(dirichlet.times, dirichlet.errorV); hold on;
% semilogy(neumann.times, neumann.errorV);
% legend('dirichlet', 'neumann')
% xlabel('time t')
% ylabel('|H|_\infty')
% 
% axes('Position',[.45 .4 .4 .3])
% box on; grid on;
% semilogy(dirichlet.times, dirichlet.errorV); hold on;
% semilogy(neumann.times, neumann.errorV);
% xlim([8.6 8.75])
% ylim([1e-4 1.05e-4])

% % fig:pde-step-time, fig:pde-step-iteration
% step8 = load('plot-data/simple-step-8.mat');
% step16 = load('plot-data/simple-step-16.mat');
% step32 = load('plot-data/simple-step-32.mat');
% figure(1);
% semilogy(step8.times, step8.errorV); hold on;
% semilogy(step16.times, step16.errorV);
% semilogy(step32.times, step32.errorV);
% xlabel('time');
% ylabel('|H|_\infty');
% title(legend('1/8, 4e-3', '1/16, 1e-3', '1/32, 4e-4'), '(h,\tau)')
% grid on;
% 
% figure(2);
% semilogy((1:numel(step8.errorV))*50, step8.errorV); hold on;
% semilogy((1:numel(step16.errorV))*50, step16.errorV);
% semilogy((1:numel(step32.errorV))*50, step32.errorV);
% xlabel('update iteration');
% ylabel('|H|_\infty');
% title(legend('1/8, 4e-3', '1/16, 1e-3', '1/32, 4e-4'), '(h,\tau)')
% grid on;

% % fig:pde-scheme-time, fig:pde-scheme-iteration
% ee = load('plot-data/explicit-euler.mat');
% ie = load('plot-data/implicit-euler.mat');
% cn = load('plot-data/crank-nicolson.mat');
% figure(1);
% semilogy(ee.times, ee.errorV); hold on;
% semilogy(ie.times, ie.errorV);
% semilogy(cn.times, cn.errorV);
% xlabel('time');
% ylabel('|H|_\infty');
% legend('explicit euler', 'implicit euler', 'crank-nicolson')
% grid on;
% 
% figure(2);
% % loglog((1:numel(ee.errorV))*50, ee.errorV); hold on;
% semilogy((1:numel(ie.errorV))*20, ie.errorV); hold on;
% semilogy((1:numel(cn.errorV))*20, cn.errorV);
% plot([100 100], ylim, 'g-.')
% xlabel('update iteration');
% ylabel('|H|_\infty');
% legend('implicit euler', 'crank-nicolson')
% grid on;

% fig:pde-adaptive-iteration
adaptive = load('plot-data/adaptive.mat');
nonadaptive = load('plot-data/non-adaptive.mat');
figure(1);
semilogy((1:numel(adaptive.errorV))*50, adaptive.errorV); hold on;
semilogy((1:numel(nonadaptive.errorV))*50, nonadaptive.errorV);
title(legend('arctan(x)', 'x'), 'space coordinate')
grid on;
xlabel('update iteration');
ylabel('|H|_\infty');

figure(2);
s = subplot(1, 2, 1); hold(s, 'on')
q1s = 0.05:0.05:0.95;
q2s = 1./sqrt(3*(1./q1s.^2-1)+1);
plot(q1s, q2s, 'r-.')
scatter(nonadaptive.sol(1, :), nonadaptive.sol(2, :), 'b.')
title('x\rightarrow x')
xlim([0 1]); ylim([0 1]);
xlabel('q_1'); ylabel('q_2')
s = subplot(1, 2, 2); hold(s, 'on')
q1s = 0.05:0.05:0.95;
q2s = 1./sqrt(3*(1./q1s.^2-1)+1);
plot(q1s, q2s, 'r-.')
scatter(adaptive.sol(1, :), adaptive.sol(2, :), 'b.')
title('x\rightarrow arctan(x)')
xlim([0 1]); ylim([0 1]);
xlabel('q_1'); ylabel('q_2')
set(gcf, 'Position', [200, 200, 560, 200])

figure(3);
s = subplot(1, 2, 1); hold(s, 'on')
plot(-2*hFunc(q1s), -2*hFunc(q2s), 'r-.')
scatter(nonadaptive.p(1, :), nonadaptive.p(2, :), 'b.')
title('x\rightarrow x')
xlim([-0.8 0]); ylim([-0.8 0]);
xlabel('p_1'); ylabel('p_2')
s = subplot(1, 2, 2); hold(s, 'on')
plot(-2*hFunc(q1s), -2*hFunc(q2s), 'r-.')
scatter(adaptive.p(1, :), adaptive.p(2, :), 'b.')
title('x\rightarrow arctan(x)')
xlim([-0.8 0]); ylim([-0.8 0]);
xlabel('p_1'); ylabel('p_2')
set(gcf, 'Position', [200, 200, 560, 200])

figure(4); hold on;
adaarc = arcsum(adaptive);
nonarc = arcsum(nonadaptive);
plot(0:(1/(numel(adaarc)-1)):1, adaarc, 'r-')
plot(0:(1/(numel(nonarc)-1)):1, nonarc, 'b-')
plot([0 1], [0 1], 'k:')
title(legend('arctan(x)', 'x', 'Location', 'Best'), 'space coordinate')
xlabel('normalized coordinate')
ylabel('normalized arc-length')
set(gcf, 'Position', [200, 200, 380, 200])

function hv = hFunc(q)
    hv = q-q.^3;
end

function arc = arcsum(s)
    dq = s.sol(:, 1:end-1) - s.sol(:, 2:end);
    dp = s.p(:, 1:end-1) - s.p(:, 2:end);
    d = sqrt(sum([dq.^2; dp.^2], 1));
    arc = cumsum(d);
    arc = arc / arc(end);
end