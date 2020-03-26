function P = Optimal_Transportcvx(X_f,Y,cost,Aeq,Aeq_1,weight)
M = size(X_f,2);
D = cost(X_f);
W = weight(X_f,Y)';
W = W./sum(W);
D1 = reshape(D,[M*M,1]);

beq = (1/M)*ones(M,1);
A = [Aeq;Aeq_1];
options = optimoptions('linprog','Algorithm','dual-simplex','display', 'off');
cvx_begin quiet
  variable x(M*M,1);
  minimize(D1'*x);
  subject to
    A*x == [beq;W'];
    x >= 0;
cvx_end

P = M.*x;
P = reshape(P,[M,M]);
cvx_status
end