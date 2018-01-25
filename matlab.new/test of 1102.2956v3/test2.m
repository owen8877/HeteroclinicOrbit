clc
clear

v0 = [0; 0; 0; 2*pi];
eps = 0.01;

lHess = notOrdinaryHessian(@bla, v0)
[V, D] = eig(lHess)
vv = V(:, 1) * eps + v0
dbla(vv)

function r = bla(v)
    a = v(1); b = v(2); c = v(3); d = v(4);
    r = sin(a)*sin(b) + cos(c+d);
end

function dblaa = dbla(v)
    dblaa = zeros(4, 1);
    eps = 1e-5;
    for i = 1:4
        hhh = zeros(4, 1); hhh(i) = eps;
        dblaa(i) = (bla(v+hhh) - bla(v-hhh)) / (2*eps);
    end
end