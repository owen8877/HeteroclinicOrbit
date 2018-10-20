function [lVal, rVal, Hqf, Hpf, Hf] = probTestGen()    
    lVal = [-1; 0; 0; 0];
    rVal = [-0.6; 0; 0; 0];
    Hqf = @Hqfunc;
    Hpf = @Hpfunc;
    Hf = @Hfunc;
end

function Hpv = Hpfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    Hpv = [q2; p2/2+df(q1)*q2];
end

function Hqv = Hqfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    Hqv = [p2*q2*ddf(q1)+2*(q2-f(q1))*df(q1); p1+p2*df(q1)-2*(q2-f(q1))];
end

function H = Hfunc(q, p)
    p1 = p(1); p2 = p(2);
    q1 = q(1); q2 = q(2);
    H = p1*q2 + p2*(p2/4+df(q1)*q2) - (q2-f(q1))^2;
end

function fv = f(x)
%     fv = 5*x^5 + 8*x^4 + 3*x^3;
    fv = -(5*x^4 + 8*x^3 + 3*x^2);
%     fv = 12*x^3 + 12*x^2;
%     fv = 2*x^3 - x/2;
%     fv = 4*x^3 - 2*x;
end

function dfv = df(x)
%     dfv = 25*x^4 + 32*x^3 + 9*x^2;
    dfv = -(20*x^3 + 24*x^2 + 6*x);
%     dfv = 36*x^2 + 24*x;
%     dfv = 6*x^2 - 1/2;
%     dfv = 12*x^2 - 2;
end

function ddfv = ddf(x)
%     ddfv = 100*x^3 + 96*x^2 + 18*x;
    ddfv = -(60*x^2 + 48*x + 6);
%     ddfv = 72*x + 24;
%     ddfv = 12*x;
%     ddfv = 24*x;
end