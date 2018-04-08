clc; clear

load data/solution.mat

% for self-made cases
lVal = [-1; -1; 0; 0];
rVal = [0; 0; 0; 0];

lAnchor = solution(:, 1);
rAnchor = solution(:, end);

resolution = 3000;
tSpan = [0 (size(solution, 2)-1)/resolution];
ntSpan = fliplr(tSpan);
fResolution = size(solution, 2)-1;

richPlotHelper(0, solution, lVal, rVal, lVal, @Hfunc, @dH, struct());
pSolution = simpleSymplecticODESolver(@qdHnv, ntSpan, rAnchor(3:4), ...
    solution(1:2, :), fResolution, struct());
solution(3:4, :) = pSolution;
richPlotHelper(0, solution, lVal, rVal, lVal, @Hfunc, @dH, struct());

%%
function vv = pdHnv(t, p, q)
    v = dHn(q, p);
    vv = v(3:4);
end

function vv = qdHnv(t, q, p)
    v = dHn(q, p);
    vv = v(1:2);
end

%%
function nv = normalize(v)
    nv = v / norm(v);
end

function dHv = dH(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
end

function dHv = nablaH(q, p)
    dHv = [Hqfunc(q, p); Hpfunc(q, p)];
end

function dHnv = dHn(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
    dHnv = dHv / norm(dHv);
end

function dHdpnv = dHdpn(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
    dHnv = dHv / norm(dHv);
    dHdpnv = dHnv(3:4);
end

function dHdqnv = dHdqn(q, p)
    dHv = [Hpfunc(q, p); -Hqfunc(q, p)];
    dHnv = dHv / norm(dHv);
    dHdqnv = dHnv(1:2);
end

%% Self-made case

function Hpv = Hpfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    Hpv = [p1+q1-q1^3; p2+q2-q2^3];
end

function Hqv = Hqfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    L = 5;
    Hqv = [p1*(1-3*q1^2)+2*((q2^2-1)*q1^2-L*(q1^2-1)*q2^2)*(2*q1*(q2^2-1)-L*2*q1*q2^2); ...
        p2*(1-3*q2^2)+2*((q2^2-1)*q1^2-L*(q1^2-1)*q2^2)*(2*q2*q1^2-L*2*q2*(q1^2-1))];
end

function H = Hfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    L = 5;
    H = (p1^2+p2^2)/2 + p1*(q1-q1^3) + p2*(q2-q2^3) + ...
        ((q2^2-1)*q1^2-L*(q1^2-1)*q2^2)^2;
end