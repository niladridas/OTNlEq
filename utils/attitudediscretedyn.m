function [enext,betanext] = attitudediscretedyn(k,e,beta,Q,T,u)
    A = Amatrix(k,T,u);
    enext = A*e;
    wbeta = mvnrnd(zeros(1,3),Q,1);
    betanext = eye(3)*beta+wbeta';
end

function A = Amatrix(k,T,U)
    u = U((k-1)*T);
    s = (T/2)*norm(u);
    p = u(1); q = u(2); r = u(3);
    Omega = [0  r -q p;
            -r  0  p q;
             q -p  0 r;
            -p -q -r 0];
    A = cos(s)*eye(4) - 0.5*T*sin(s)*Omega/s;
end