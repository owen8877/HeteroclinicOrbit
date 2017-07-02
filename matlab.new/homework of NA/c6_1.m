depthSpan = 1:8;

xRange = [1 2];
y0 = 0.4;

f = @forwardf;
bf = @backwardf;

semiSolution = ode23(f, xRange, y0);
% h = 0.1;
% for exponent = 1:7
%     ooo = ode45(f, xRange, y0+2^(-exponent+1)*h);
%     fprintf('%d\t%e\n', exponent, deval(ode45Solution, xRange(2)) - deval(ooo, xRange(2)));
% end
disp('ode23')
preciseY = deval(semiSolution, xRange(2));
fprintf('%e\t\t%f\t%f\n', nan, deval(semiSolution, xRange(2)), deval(semiSolution, xRange(2))- preciseY);

ODETester('naive forward method', @odeEulerForwardMethod, depthSpan, f, xRange, y0);
ODETester('naive backward method', @odeEulerBackwardMethod, depthSpan, bf, xRange, y0);
ODETester('Kutta 3 method', @ode3KuttaMethod, depthSpan, f, xRange, y0);
ODETester('Heun 3 method', @ode3HeunMethod, depthSpan, f, xRange, y0);
ODETester('Kutta 4 method', @ode4KuttaMethod, depthSpan, f, xRange, y0);
ODETester('LinearMulti 3 Explicit method', @ode3LinearMultiMethod, depthSpan, f, xRange, y0);
ODETester('Improved Euler Method', @odeImprovedEulerMethod, depthSpan, f, xRange, y0);
ODETester('Predict Milne 3 Method', @ode3MilneMethod, depthSpan, f, xRange, y0);
ODETester('Predict Adams 4 Method', @ode4AdamsMethod, depthSpan, f, xRange, y0);

function yd = forwardf(x, y, ~)
    yd = -(1/x*x) - y/x - y*y;
end

function yd = backwardf(x, y, h)
    yd = (-(1/x*x) - y/x - y*y) * (1 - h * (2*y+1/x) / 2);
end

function yd = exampleff(x, y, ~)
    yd = x^3 - y/x;
end

function yd = examplebf(x, y, ~)
    yd = x^3 - y/x;
end

function ODETester(desc, ODEMethod, depthSpan, f, xRange, y0)
    h = 0.1;
    disp(desc)
    resulrArray = zeros(size(depthSpan, 1));
    for exponent = depthSpan
        solution = ODEMethod(f, xRange, y0, (xRange(2) - xRange(1))/(2^(-exponent+1) * h)+1);
        yf = solution(end, 2);
        % fprintf('%e\t%f\t%e\n', 2^(-exponent+1) * h, yf, yf-preciseY);
        resulrArray(exponent-depthSpan(1)+1) = yf;
    end
    rErrorTermArray = resulrArray(2:end) - resulrArray(1:end-1);
    coeff = polyfit(depthSpan(2:end), log2(rErrorTermArray), 1);
    fprintf('order : %f\n\n', -coeff(1));
end

function solutionArray = ODESolverWrapper(core, f, xRange, y0, steps)
    x0 = xRange(1);
    xf = xRange(2);
    xs = linspace(x0, xf, steps);
    dx = (xf - x0) / (steps-1);
    ys = zeros(steps-1, 1);
    ys(1) = y0;
    ys = core(steps, dx, xs, ys, f);
    solutionArray = [xs' ys];
end

function solutionArray = odeEulerForwardMethod(f, xRange, y0, steps)
    function ys = core(steps, dx, xs, ys, f)
        for i = 1:(steps-1)
            x = xs(i);
            ys(i+1) = ys(i) + dx * f(x, ys(i), dx);
        end
    end

    solutionArray = ODESolverWrapper(@core, f, xRange, y0, steps);
end

function solutionArray = odeEulerBackwardMethod(f, xRange, y0, steps)
    function ys = core(steps, dx, xs, ys, f)
        for i = 1:(steps-1)
            x = xs(i+1);
            ys(i+1) = ys(i) + dx * f(x, ys(i), dx);
        end
    end

    solutionArray = ODESolverWrapper(@core, f, xRange, y0, steps);
end

function solutionArray = ode3KuttaMethod(f, xRange, y0, steps)
    function ys = core(steps, dx, xs, ys, f)
        for i = 1:(steps-1)
            x = xs(i);
            y = ys(i);
            K1 = f(x, y);
            K2 = f(x+dx/2, y+dx*K1/2);
            K3 = f(x+dx, y-dx*(K1-2*K2));
            ys(i+1) = ys(i) + dx/6 * (K1 + K2*4+K3);
        end
    end

    solutionArray = ODESolverWrapper(@core, f, xRange, y0, steps);
end

function solutionArray = ode3HeunMethod(f, xRange, y0, steps)
    function ys = core(steps, dx, xs, ys, f)
        for i = 1:(steps-1)
            x = xs(i);
            y = ys(i);
            K1 = f(x, y);
            K2 = f(x+dx/3, y+dx*K1/3);
            K3 = f(x+dx*2/3, y+dx*2*K2/3);
            ys(i+1) = ys(i) + dx/4 * (K1 + K3*3);
        end
    end

    solutionArray = ODESolverWrapper(@core, f, xRange, y0, steps);
end

function solutionArray = ode4KuttaMethod(f, xRange, y0, steps)
    function ys = core(steps, dx, xs, ys, f)
        for i = 1:(steps-1)
            x = xs(i);
            y = ys(i);
            K1 = f(x, y);
            K2 = f(x+dx/3, y+dx*K1/3);
            K3 = f(x+dx*2/3, y+dx*(-K1/3+K2));
            K4 = f(x+dx, y+dx*(K1-K2+K3));
            ys(i+1) = ys(i) + dx/8 * (K1 + K2*3 + K3*3 + K4);
        end
    end

    solutionArray = ODESolverWrapper(@core, f, xRange, y0, steps);
end

function solutionArray = odeImprovedEulerMethod(f, xRange, y0, steps)
    function ys = core(steps, dx, xs, ys, f)
        for i = 1:(steps-1)
            x = xs(i);
            y = ys(i);
            y_ = y + dx * f(x, y);
            ys(i+1) = y + dx/2 * (f(x, y) + f(x+dx, y_));
        end
    end

    solutionArray = ODESolverWrapper(@core, f, xRange, y0, steps);
end

function solutionArray = ode3LinearMultiMethod(f, xRange, y0, steps)
    function ys = core(steps, dx, xs, ys, f)
        ys(2) = ys(1) + dx * f(xs(1), ys(1));
        ys(3) = ys(2) + dx / 2 * (3*f(xs(2), ys(2)) - f(xs(1), ys(1)));

        for i = 3:(steps-1)
            P1 = f(xs(i), ys(i));
            P2 = f(xs(i-1), ys(i-1));
            P3 = f(xs(i-2), ys(i-2));
            ys(i+1) = ys(i-1) + dx/3 * (7*P1 - 2*P2 + P3);
        end
    end

    solutionArray = ODESolverWrapper(@core, f, xRange, y0, steps);
end

function solutionArray = ode3MilneMethod(f, xRange, y0, steps)
    function ys = core(steps, dx, xs, ys, f)
%         ys(2) = ys(1) + dx * f(xs(1), ys(1));
%         ys(3) = ys(2) + dx/2 * (3*f(xs(2), ys(2)) - f(xs(1), ys(1)));
%         ys(4) = ys(2) + dx/3 * (7*f(xs(3), ys(3) - 2*f(xs(2), ys(2)) + f(xs(1), ys(1))));
        for i = 1:3
            x = xs(i);
            y = ys(i);
            K1 = f(x, y);
            K2 = f(x+dx/3, y+dx*K1/3);
            K3 = f(x+dx*2/3, y+dx*(-K1/3+K2));
            K4 = f(x+dx, y+dx*(K1-K2+K3));
            ys(i+1) = ys(i) + dx/8 * (K1 + K2*3 + K3*3 + K4);
        end

        for i = 4:(steps-1)
            y_ = ys(i-3) + dx/3 * (8*f(xs(i), ys(i)) - 4*f(xs(i-1), ys(i-1)) + 8*f(xs(i-2), ys(i-2)));
            P1 = f(xs(i+1), y_);
            P2 = f(xs(i), ys(i));
            P3 = f(xs(i-1), ys(i-1));
            ys(i+1) = ys(i-1) + dx/3 * (P1 + 4*P2 + P3);
        end
    end

    solutionArray = ODESolverWrapper(@core, f, xRange, y0, steps);
end

function solutionArray = ode4AdamsMethod(f, xRange, y0, steps)
    function ys = core(steps, dx, xs, ys, f)
%         ys(2) = ys(1) + dx * f(xs(1), ys(1));
%         ys(3) = ys(2) + dx/2 * (3*f(xs(2), ys(2)) - f(xs(1), ys(1)));
%         ys(4) = ys(2) + dx/3 * (7*f(xs(3), ys(3) - 2*f(xs(2), ys(2)) + f(xs(1), ys(1))));
        for i = 1:3
            x = xs(i);
            y = ys(i);
            K1 = f(x, y);
            K2 = f(x+dx/3, y+dx*K1/3);
            K3 = f(x+dx*2/3, y+dx*(-K1/3+K2));
            K4 = f(x+dx, y+dx*(K1-K2+K3));
            ys(i+1) = ys(i) + dx/8 * (K1 + K2*3 + K3*3 + K4);
        end

        for i = 4:(steps-1)
            P0 = f(xs(i), ys(i));
            P1 = f(xs(i-1), ys(i-1));
            P2 = f(xs(i-2), ys(i-2));
            P3 = f(xs(i-3), ys(i-3));
            y_ = ys(i) + dx/24 * (55*P0 - 59*P1 + 37*P2 - 9*P3);
            ys(i+1) = ys(i) + dx/24 * (9*f(xs(i+1), y_) + 19*P0 - 5*P1 + P2);
        end
    end

    solutionArray = ODESolverWrapper(@core, f, xRange, y0, steps);
end