function fltrd_data =  OT_filtertrial1(theta_samples,measured_output,cost,weight,OT_constants,Optimal_Transport,L)
    M = size(theta_samples,2);
    [Aeq,Aeq_1] = OT_constants(M);
    P = Optimal_Transport(theta_samples,measured_output,cost,Aeq,Aeq_1,weight);
%     fltrd_data = theta_samples*P;
    % The P matrix is perturbed to come closer to the equality constraint
    % for the pendulum problem
    X1 = theta_samples(1:2,:);
    Q = X1'*X1;
    A = P'*Q*P;
    B = 2*P'*Q;
    C = L^2*rand(size(P,1),size(P,1));
    c1 = diag(C);
    c1 = L^2-c1;
    C = C + diag(c1);
    delP = pinv(B)*(C - A);
    P = P + delP;
    fltrd_data = theta_samples*P;
end