function richPlotHelper(itr, solution, lVal, rVal, rVal2, Hfunc, dH, options)
    if ~isfield(options, 'extraPlot')
        options.extraPlot = false;
    end
    if ~isfield(options, 'L')
        options.L = 5;
    end
    
    figure
    plotHelper(subplot(2, 2, 1), solution([1 3], :), 'q_1', 'p_1', [lVal([1 3]) rVal([1 3]) rVal2([1 3], :)]);
    plotHelper(subplot(2, 2, 2), solution([2 4], :), 'q_2', 'p_2', [lVal([2 4]) rVal([2 4]) rVal2([2 4], :)]);
    plotHelper(subplot(2, 2, 3), solution([1 2], :), 'q_1', 'q_2', [lVal([1 2]) rVal([1 2]) rVal2([1 2], :)]);
    plotHelper(subplot(2, 2, 4), solution([3 4], :), 'p_1', 'p_2', [lVal([3 4]) rVal([3 4]) rVal2([3 4], :)]);
    
    if options.extraPlot
        q2 = linspace(-0.99, -0.01, 99);
        K = options.L * q2.^2 ./ (1 - q2.^2);
        q1 = -sqrt(K./(1+K));
        p1 = 2*(q1.^3-q1);
        p2 = 2*(q2.^3-q2);
        plotHelper(subplot(2, 2, 1), [q1; p1], 'q_1', 'p_1', []);
        plotHelper(subplot(2, 2, 2), [q2; p2], 'q_2', 'p_2', []);
        plotHelper(subplot(2, 2, 3), [q1; q2], 'q_1', 'q_2', []);
        plotHelper(subplot(2, 2, 4), [p1; p2], 'p_1', 'p_2', []);
    end
%     title(sprintf('Itr=%d', itr))

    H = zeros(1, size(solution, 2));
    for i = 1:size(solution, 2)
        H(i) = Hfunc(solution(1:2, i), solution(3:4, i));
    end

    figure
    hold on; grid on
    plot(H)
    line(get(gca,'XLim'), [0 0], 'Color', 'r');
    
%     title(sprintf('H Itr=%d', itr))
    xlabel('node')
    ylabel('|H| error')
    
    if options.extraPlot
        q2 = solution(2, :);
        K = options.L * q2.^2 ./ (1 - q2.^2);
        q1 = -sqrt(K./(1+K));
        p1 = 2*(q1.^3-q1);
        p2 = 2*(q2.^3-q2);
        diff = solution([1 3 4], :) - [q1; p1; p2];
        
        figure
        plot(sqrt(sum(diff.^2, 1)))
        title('Estimated error in trajectory.')
    end
% 
%     errorIndH = zeros(1, size(solution, 2)-1);
%     for i = 1:size(solution, 2)-1
%         dqp = solution(:, i+1) - solution(:, i);
%         dHdqp = dH(solution(1:2, i), solution(3:4, i));
%         dqp = dqp / norm(dqp) * norm(dHdqp);
%         errorIndH(i) = norm(dqp - dHdqp);
%     end
% 
%     figure
%     plot(errorIndH)
%     title(sprintf('errorIndH Itr=%d', itr))
end

