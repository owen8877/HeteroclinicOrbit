function C = solutionLength(solution)
    C = 0;
    for i = 1:size(solution, 2)-1
        C = C + norm(solution(:, i) - solution(:, i+1));
    end
end