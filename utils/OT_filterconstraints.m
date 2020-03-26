function fltrd_data =  OT_filterconstraints(theta_samples,measured_output,cost,weight,OT_constants,Optimal_Transport,L)
    M = size(theta_samples,2);
    [Aeq,Aeq_1] = OT_constants(M);
    P = Optimal_Transport(theta_samples,measured_output,cost,Aeq,Aeq_1,weight,L);
    fltrd_data = theta_samples*P;
end