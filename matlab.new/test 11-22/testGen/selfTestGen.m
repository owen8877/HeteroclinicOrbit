function [lVal, rVal, Hqf, Hpf, Hf] = selfTestGen()
    lVal = [-1; -1; 0; 0];
    rVal = [0; 0; 0; 0];
    Hqf = @Hqfunc;
    Hpf = @Hpfunc;
    Hf = @Hfunc;
end

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