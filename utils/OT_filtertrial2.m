function fltrd_data =  OT_filtertrial2(theta_samples,measured_output,cost,weight,OT_constants,Optimal_Transport,L)
    M = size(theta_samples,2);
    [Aeq,Aeq_1] = OT_constants(M);
    P = Optimal_Transport(theta_samples,measured_output,cost,Aeq,Aeq_1,weight);
%     fltrd_data = theta_samples*P;
    % The P matrix is perturbed to come closer to the equality constraint
    % for the pendulum problem
    X1 = theta_samples(1:2,:);
    Q = X1'*X1;
    Qhalf = sqrtm(Q);
    np = size(P,1);
    A = L^2*eye(np);
    A1 = A - P'*Q*P;
    B = 2*P'*Q;
    cvx_begin sdp
     variable delP(np,np) %symmetric
      minimize(norm(delP,'fro'))
      subject to
     [A1-B*delP delP'*Qhalf'; Qhalf*delP eye(np)] >= 0;
    cvx_end
    P = P + delP;
    fltrd_data = theta_samples*P;
end