function P = Optimal_Transportcvxconstraints(X_f,Y,cost,Aeq,Aeq_1,weight,L)
M   = size(X_f,2);
D   = cost(X_f);
W   = weight(X_f,Y)';
W   = W./sum(W);
D1  = reshape(D,[M*M,1]);
beq = (1/M)*ones(M,1);
A   = [Aeq;Aeq_1];
Q   = X_f'*X_f;
cvx_begin sdp
    variable x(M*M,1);
    minimize(D1'*x);
    subject to
        A*x == [beq;W'];
        x >= 0;
        for i = 1:M
            Mtmp = zeros(M,M*M);
            Mtmp(:,M*(i-1)+1:M*i) = eye(M,M);
            Qtmp = Mtmp'*Q*Mtmp;
            Qtmp = 0.5*(Qtmp + Qtmp');
            X1 = [L^2/(M*M), x'*sqrtm(Qtmp);sqrtm(Qtmp)*x ,eye(M*M,M*M)];
            X1 >= 0;
        end   
cvx_end

P = M.*x;
P = reshape(P,[M,M]);
cvx_status
end