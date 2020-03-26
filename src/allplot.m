clc;clear;close;
L        = 1.0;
%% Plotting the RMSE error of the estimates
load('../data/Datapendulum.mat');

% Estimate across all the MC run and across time
x_only    = XA_meanonlyOT;
x_proj    = XA_meanOTProj;
x_ma      = XA_meanOTMA;
x_nleq    = XA_meanOTNLEQ;
x_nleqma  = XA_meanOTNLEQMA;

% Cov across all the MC run and across time
xc_only    = XA_covonlyOT;
xc_proj    = XA_covOTProj;
xc_ma      = XA_covOTMA;
xc_nleq    = XA_meanOTNLEQ;
xc_nleqma  = XA_covOTNLEQM;

% Estimate of x and y across all the MC run and across time
xy_only     =  XA_meanonlyOT(:,1:2,:);
xy_proj     = XA_meanOTProj(:,1:2,:);
xy_ma       = XA_meanOTMA(:,1:2,:);
xy_nleq     = XA_meanOTNLEQ(:,1:2,:);
xy_nleqma   = XA_meanOTNLEQMA(:,1:2,:);

% Cov of x and y across all the MC run and across time
xyc_only     = XA_covonlyOT(:,1:2,:);
xyc_proj     = XA_covOTProj(:,1:2,:);
xyc_ma       = XA_covOTMA(:,1:2,:);
xyc_nleq     = XA_covOTNLEQ(:,1:2,:);
xyc_nleqma   = XA_covOTNLEQM(:,1:2,:);

% Average of x and y over MC runs
xycsampavg_only     = mean(xyc_only,3);
xycsampavg_proj     = mean(xyc_proj,3);
xycsampavg_ma       = mean(xyc_ma,3);
xycsampavg_nleq     = mean(xyc_nleq,3);
xycsampavg_nleqma   = mean(xyc_nleqma,3);

% Average of x and y first over MC runs and then over time
xyctimeavg_only     = mean(xycsampavg_only,1);
xyctimeavg_proj     = mean(xycsampavg_proj,1);
xyctimeavg_ma       = mean(xycsampavg_ma,1);
xyctimeavg_nleq     = mean(xycsampavg_nleq,1);
xyctimeavg_nleqma   = mean(xycsampavg_nleqma,1);



% Constraint value across MC runs with time 
L_only     =  sqrt(squeeze(sum(xy_only.*xy_only,2)));
L_proj     = sqrt(squeeze(sum(xy_proj.*xy_proj,2)));
L_ma       = sqrt(squeeze(sum(xy_ma.*xy_ma,2)));
L_nleq     = sqrt(squeeze(sum(xy_nleq.*xy_nleq,2)));
L_nleqma   = sqrt(squeeze(sum(xy_nleqma.*xy_nleqma,2)));

nt = size(L_only,1); ns = size(L_only,2);

% Constraint error across MC runs with time 
E_only     = L_only - L*ones(nt,ns);
E_proj     = L_proj - L*ones(nt,ns);
E_ma       = L_ma - L*ones(nt,ns);
E_nleq     = L_nleq - L*ones(nt,ns);
E_nleqma   = L_nleqma - L*ones(nt,ns);


% RMSE error accross MC runs with time
rmse_only    = sqrt((1/ns)*sum(E_only.*E_only,2));
rmse_proj    = sqrt((1/ns)*sum(E_proj.*E_proj,2));
rmse_ma      = sqrt((1/ns)*sum(E_ma.*E_ma,2));
rmse_nleq    = sqrt((1/ns)*sum(E_nleq.*E_nleq,2));
rmse_nleqma  = sqrt((1/ns)*sum(E_nleqma.*E_nleqma,2));


% Average RMSE error
avgrmse_only    = sum(rmse_only)/nt;
avgrmse_proj    = sum(rmse_proj)/nt;
avgrmse_ma      = sum(rmse_ma)/nt;
avgrmse_nleq    = sum(rmse_nleq)/nt;
avgrmse_nleqma  = sum(rmse_nleqma)/nt;

% Individual RMSE error

% RMSE for position

% RMSE for velocity


% Trace of error covariance matrix for position


% plot(T,log10(rmse_only)); hold on;
% plot(T,log10(rmse_proj)); hold on;
% plot(T,log10(rmse_ma)); hold on;
% plot(T,log10(rmse_nleq)); hold on;
% plot(T,log10(rmse_nleqma)); 

Tp = 2*pi*sqrt(L/9.8);
% 
semilogy(T/Tp,(rmse_only)); hold on;
semilogy(T/Tp,(rmse_proj)); hold on;
semilogy(T/Tp,(rmse_ma)); hold on;
semilogy(T/Tp,(rmse_nleq)); hold on;
semilogy(T/Tp,(rmse_nleqma)); 
legend('OT','OTProj','OTMA','OTNLeq','OTNLeqMA');