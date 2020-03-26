function de = attitudedyn(t,e,u)
%ATTITUDESYN Summary of this function goes here
%   Detailed explanation goes here
% Attitude is defined using a quaternion
% u in the angular velocity used as an input
u = u(t);
p = u(1); q = u(2); r = u(3);
Omega = [0  r -q p;
        -r  0  p q;
         q -p  0 r;
        -p -q -r 0];
     
de = -0.5*Omega*e;
end

