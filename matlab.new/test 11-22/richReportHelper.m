function richReportHelper(itr, solution, lVal, rVal, rVal2, Hfunc, dH)
    H = zeros(1, size(solution, 2));
    for i = 1:size(solution, 2)
        H(i) = Hfunc(solution(1:2, i), solution(3:4, i));
    end
    fprintf('H err\t%.3e\n', max(abs(H)))

    errorIndH = zeros(1, size(solution, 2)-1);
    for i = 1:size(solution, 2)-1
        dqp = solution(:, i+1) - solution(:, i);
        dHdqp = dH(solution(1:2, i), solution(3:4, i));
        dqp = dqp / norm(dqp) * norm(dHdqp);
        errorIndH(i) = norm(dqp - dHdqp);
    end
    fprintf('dH err\t%.3e\n', max(abs(errorIndH)))
end

