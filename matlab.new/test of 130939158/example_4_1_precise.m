B = [[-1 4]; [-4 -1]];
% B = [[-1 0.5]; [-0.5 -1]];

t = linspace(0, 50, 10000);

pPrecise = zeros(2, size(t, 2));
xPrecise = zeros(2, size(t, 2));
p0 = [4; 1];
for i = 1:size(t, 2)
    pPrecise(:, i) = expm(B' * t(i)) * p0;
    xPrecise(:, i) = (B + B') \ pPrecise(:, i);
end

figure
subplot(1, 2, 1);
plot(pPrecise(1, :), pPrecise(2, :));
subplot(1, 2, 2);
plot(xPrecise(1, :), xPrecise(2, :));