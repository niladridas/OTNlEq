clc;clear;close;
L        = 1.0;
%% Plotting the RMSE error of the estimates
load('../data/Xreal');
load('../data/Datapendulum1.mat');

%% RMSE Error
Xall = repmat(X,1,1,50);
EonlyOT    = XA_meanonlyOT - Xall; E2onlyOT  = EonlyOT.*EonlyOT;
EOTMA      = XA_meanOTMA - Xall;   E2OTMA    = EOTMA.*EOTMA; 
EOTNLEQ    = XA_meanOTNLEQ - Xall; E2OTNLEQ  = EOTNLEQ.*EOTNLEQ;
EOTProj    = XA_meanOTProj - Xall; E2OTProj  = EOTProj.*EOTProj;
EOTNLEQMA  = XA_meanOTNLEQMA - Xall; E2OTNLEQMA = EOTNLEQMA.*EOTNLEQMA;
E2onlyOT = sqrt(mean(sum(E2onlyOT(:,1:2,:),2),3));
E2OTMA   = sqrt(mean(sum(E2OTMA(:,1:2,:),2),3));
E2OTNLEQ = sqrt(mean(sum(E2OTNLEQ(:,1:2,:),2),3));
E2OTProj = sqrt(mean(sum(E2OTProj(:,1:2,:),2),3));
E2OTNLEQMA = sqrt(mean(sum(E2OTNLEQMA(:,1,:),2),3));

% RMSE along Monte Carlo run
plot(T,E2onlyOT); hold on;
plot(T,E2OTMA);hold on;
plot(T,E2OTNLEQ); hold on;
plot(T,E2OTProj);hold on;
plot(T,E2OTNLEQMA);
legend('OTF','OTMA','OTNLeq','OTProj','OTNLEQMA')
%% Average RMSE Error


% % Average across MC run
% XonlyOT    = mean(XA_meanonlyOT(:,1:2,:),3) ;
% XOTMA      = mean(XA_meanOTMA(:,1:2,:),3) ;
% XOTNLEQ    = mean(XA_meanOTNLEQ(:,1:2,:),3) ;
% XOTProj    = mean(XA_meanOTProj(:,1:2,:),3) ;
% XOTNLEQMA  = mean(XA_meanOTNLEQMA(:,1:2,:),3);
% 
% 
% plot(T,XOTNLEQMA(:,2)); hold on; plot(T,Xall(:,2))
