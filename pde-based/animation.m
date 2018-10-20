clear; %clc
addpath util
global counter

AB = load('backup/matlabAB.mat');
BC = load('backup/matlabBC.mat');
CD = load('backup/matlabCD.mat');
AD = load('backup/matlabAD.mat');

figure(1)

counter = 0;
% display(AB, 0.01)
% display(BC, 0.01)
% display(CD, 0.01)
display(AD, 0.01)
% AB.sol = rotate(AB.sol);
% display(AB, 0.01)

% test.sol = zeros(6, 100);
% for i = 1:100
%     for j = 1:3
%         test.sol(:, i) = [ ...
%             cos(i/50); sin(i/50); ...
%             0; 0; ...
%             -cos(i/50); -sin(i/50); ...
%         ];
%     end
% end
% test.sol = rotate(test.sol);
% display(test, 0.01)

function display(data, dt)
    global counter
    for i = 1:size(data.sol, 2)
        clf
        for j = 1:size(data.sol, 1)/2
            circle(data.sol(2*j-1, i), data.sol(2*j, i), 2^(1/6)/2);
        end
        axis([-2.5 2.5 -2.5 2.5])
        drawnow
        counter = counter + 1;
%         saveas(gcf, sprintf('output/%d.png', counter), 'png')
        pause(dt)
    end
end

function h = circle(x, y, r)
    d = r*2;
    px = x-r;
    py = y-r;
    h = rectangle('Position', [px py d d], 'Curvature', [1, 1]);
    daspect([1, 1, 1])
end