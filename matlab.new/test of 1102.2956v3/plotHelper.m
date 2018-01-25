function plotHelper(ax, data, xl, yl, interest)

xData = data(1, :);
yData = data(2, :);

hold(ax, 'on');
plot(xData, yData);
xlabel(xl);
ylabel(yl);

plot(xData(1), yData(1), 'ro');
plot(xData(end), yData(end), 'bo');

if numel(interest) > 0
    plot(interest(1, :), interest(2, :), 'k+');
end

end

