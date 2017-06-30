depthSpan = 2:6;

ode45Solution = ode45(@forwardf, [1 2], -1);
disp('ode45')
preciseY = deval(ode45Solution, 2);
fprintf('%e\t\t%f\t%f\n', nan, preciseY, 0);

disp('naive forward method')
for exponent = depthSpan
    naiveForward = odeEulerMethod(@forwardf, [1 2], -1, 10^exponent - 1, true);
    fprintf('%e\t%f\t%f\n', 10^(-exponent), naiveForward(end, 2), naiveForward(end, 2)-preciseY);
end

disp('naive backward method')
for exponent = depthSpan
    naiveBackward = odeEulerMethod(@backwardf, [1 2], -1, 10^exponent - 1, false);
    fprintf('%e\t%f\t%f\n', 10^(-exponent), naiveBackward(end, 2), naiveBackward(end, 2)-preciseY);
end

disp('naive Kutta 3 method')
for exponent = depthSpan
    naive = ode3KuttaMethod(@forwardf, [1 2], -1, 10^exponent - 1);
    fprintf('%e\t%f\t%f\n', 10^(-exponent), naive(end, 2), naive(end, 2)-preciseY);
end

function yd = forwardf(x, y, ~)
    yd = -(1/x*x) - y/x - y*y;
end

function yd = backwardf(x, y, h)
    yd = (-(1/x*x) - y/x - y*y) * (1 - h * (2*y+1/x) / 2);
end

function solutionArray = odeEulerMethod(f, xRange, y0, steps, forwardOn)
    x0 = xRange(1);
    xf = xRange(2);
    xs = linspace(x0, xf, steps);
    dx = (xf - x0) / steps;
    ys = zeros(steps, 1);
    ys(1) = y0;
    for i = 1:(steps-1)
        if forwardOn
            x = xs(i);
        else
            x = xs(i+1);
        end
        ys(i+1) = ys(i) + dx * f(x, ys(i), dx);
    end
    solutionArray = [xs' ys];
end

function solutionArray = ode3KuttaMethod(f, xRange, y0, steps)
    x0 = xRange(1);
    xf = xRange(2);
    xs = linspace(x0, xf, steps);
    dx = (xf - x0) / steps;
    ys = zeros(steps, 1);
    ys(1) = y0;
    for i = 1:(steps-1)
        x = xs(i);
        y = ys(i);
        K1 = f(x, y);
        K2 = f(x+dx/2, y+dx*K1/2);
        K3 = f(x+dx, y-dx*(K1-2*K2));
        ys(i+1) = ys(i) + dx/6 * (K1 + K2*4+K3);
    end
    solutionArray = [xs' ys];
end