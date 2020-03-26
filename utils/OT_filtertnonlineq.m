function fltrd_data =  OT_filtertnonlineq(theta_samples,measured_output,cost,weight,OT_constants,Optimal_Transport,L)
    M = size(theta_samples,2);
    [Aeq,Aeq_1] = OT_constants(M);
    P = Optimal_Transport(theta_samples,measured_output,cost,Aeq,Aeq_1,weight);
    fltrd_data = theta_samples*P;
    Ns = size(fltrd_data,2);
    Nd = size(fltrd_data,1);
    %Pendulum specific
    X1 = fltrd_data(1:2,:);
    X1mean = mean(fltrd_data,2);
    X1 = sum(X1.*X1,1);
    Edk = mean(X1);
    Sigdd = cov(X1');
    Sigxd = (1/(Ns-1))*sum((fltrd_data-X1mean.*ones(Nd,Ns)).*(X1-Edk),2);
    Kp = Sigxd/Sigdd;
    fltrd_data = fltrd_data + Kp*(L^2*ones(1,Ns)-X1);
end