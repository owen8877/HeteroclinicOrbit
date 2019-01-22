function plotHelper(ax, data, xl, yl, interest)
    xData = data(1, :);
    yData = data(2, :);

    hold(ax, 'on'); grid on
    plot(xData, yData);
    xlabel(xl);
    ylabel(yl);

    plot(xData(1), yData(1), 'ro');
    plot(xData(end), yData(end), 'bo');

    if numel(interest) >= 1
        plot(interest(1, 1), interest(2, 1), 'm^', 'MarkerFaceColor', 'm');
    end
    if numel(interest) >= 2
        plot(interest(1, 2), interest(2, 2), 'm^', 'MarkerFaceColor', 'c');
    end
    if numel(interest) >= 3
        plot(interest(1, 2), interest(2, 2), 'm^', 'MarkerFaceColor', 'k');
    end
end

