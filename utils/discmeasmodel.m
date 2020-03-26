%% Measurement Model
% Measurements are available at a much lower rate of 10*Tsample
function y = discmeasmodel(E,R)
   b = E(1); c = E(2); d = E(3); a = E(4);
   C = [2*a^2-1+2*b^2, 2*(b*c+a*d), 2*(b*d-a*c);
        2*(b*c-a*d), 2*(a^2+c^2)-1, 2*(c*d+a*b);
        2*(b*d+a*c), 2*(c*d-a*b), 2*(a^2+d^2)-1];
   v =   mvnrnd(zeros(1,6),R,1);
   y1 = C*[1;0;0];
   y2 = C*[0;1;0];
   y = [y1;y2]+v';
   y = y';
end