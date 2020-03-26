clc;close;clear;
rng('default');

Ns      = 10;      % Sample size
disp1   = 0;       % Plot flag 1
disp2   = 1;       % Plot flag 2
%
g        = 9.8; 
L        = 1.0;
tinitial = 0;
tfinal   = 5.0;
X0       = [L*cos(20*pi/180),L*sin(20*pi/180),0.0,0.0]';
deltat   = 0.05;
tspan    = tinitial:deltat:tfinal;
params.g = g;
params.L = L;
dyn_pendulum = @(t,x) pendulumfixedlengthdyn(t,x,params);
options = odeset('RelTol',1e-12,'AbsTol',1e-10);
[T,X]    = ode45(dyn_pendulum,tspan,X0,options);
% Proof of Numerical error 
% PLOT 1
if disp1 ==1
    figure(1);subplot(2,2,1);scatter(X(:,1),-1*X(:,2),'filled','r') ; xlabel('X axis'); ylabel('Y axis');    grid on;
    subplot(2,2,2); plot(T,sqrt(sum(X(:,1:2).*X(:,1:2),2))-L); ylabel('Length error in the pendulum'); xlabel('Time [secs]');    grid on;
    subplot(2,2,3); plot(T,X(:,1));  ylabel('X position'); xlabel('Time [secs]');    grid on;
    subplot(2,2,4); plot(T,X(:,2)) ; ylabel('Y position'); xlabel('Time [secs]');    grid on;
end
%%
Nmeas = size(X,1);
% Generate measurements
% Only position measurements
% Noise S.D. is 0.5 meters in each
meas_noise_cov = [0.5,0.5].^2;
% Generate random noise and add to each of the 2 channels
meas_noise_mu =  [0,0];
meas_noise    = mvnrnd(meas_noise_mu, diag(meas_noise_cov),Nmeas);
Y_meas        = X(:,1:2) + meas_noise;
%%
% Generate initial samples satisfying the constraint.
theta_0  = 30*pi/180; 
theta_0width = 20*pi/180;   % +/- 5 degress about theta_0;
theta_0 = theta_0 + theta_0width*2*(rand(1,Ns)-0.5);
X0 = [L*cos(theta_0);L*sin(theta_0);zeros(2,Ns)];
%%
% OT filtering
X_prior = X0;
cost    = @(X) distance_matrix(X);
% weight  = @(X,Y) mvnpdf(Y-X(1:2,:)',meas_noise_mu,diag(meas_noise_cov));
weight  = @(X,Y) weightPM(X,Y,meas_noise_cov,meas_noise_mu,L);

X_mean  = zeros(Nmeas,4);
X_cov   = zeros(Nmeas,4); % sig_xx, sig_yy, sig_velxx, sig_velyy
for i = 1:Nmeas
    % UPDATE
%     X_post = OT_filter(X_prior,Y_meas(i,:),cost,weight,@OT_constants,@Optimal_Transport);
    X_post = OT_filtertnonlineq(X_prior,Y_meas(i,:),cost,weight,@OT_constants,@Optimal_Transport,L);
%     [X_post,X_postP] =  OT_filtertnonlineqproj(X_prior,Y_meas(i,:),cost,weight,@OT_constants,@Optimal_Transport,L);
%       X_post = OT_filtertrial2(X_prior,Y_meas(i,:),cost,weight,@OT_constants,@Optimal_Transport,L);
%     X_post   = OT_filterconstraints(X_prior,Y_meas(i,:),cost,weight,@OT_constants,@Optimal_Transportcvxconstraints,L);
    msgfilt  = sprintf('%i th filtering step of total %i steps',i,Nmeas);
    disp(msgfilt);
    %% Only for projection 
%     X_mean(i,:) = mean(X_postP,2)';
%     X_cov(i,:)  = diag(cov(X_postP'))';
    
    X_mean(i,:) = mean(X_post,2)';
    X_cov(i,:)  = diag(cov(X_post'))';

    % Save mean and cov of filtered states
    % PROPAGATION
    if i==length(T)
        break;
    end
    tspan    = [T(i),T(i+1)];
    for j = 1:Ns
        [~,Xtmp]     =  ode45(dyn_pendulum,tspan,X_post(:,j),options);
        X_prior(:,j) =  Xtmp(end,:)';
    end
end
%% PLOT 2
if disp2 ==1
    figure(2);subplot(2,3,1);scatter(X_mean(:,1),-1*X_mean(:,2),'filled','r') ; xlabel('X axis'); ylabel('Y axis');    grid on;
    subplot(2,3,4); plot(T,abs(sqrt(sum(X_mean(:,1:2).*X_mean(:,1:2),2))-L)); ylabel('Length error in the pendulum'); xlabel('Time [secs]');    grid on;
    subplot(2,3,2); plot(T,X(:,1),T,X_mean(:,1));  ylabel('X position'); xlabel('Time [secs]'); legend('real','filtered');   grid on;
    subplot(2,3,3); plot(T,X(:,2),T,X_mean(:,2)) ; ylabel('Y position'); xlabel('Time [secs]');  legend('real','filtered');  grid on;
    subplot(2,3,5); plot(T,sqrt(X_cov(:,1))) ; ylabel('X position SD'); xlabel('Time [secs]');  grid on;
    subplot(2,3,6); plot(T,sqrt(X_cov(:,2))) ; ylabel('Y position SD'); xlabel('Time [secs]');  grid on;    
end

function W = weightPM(X_f,Ya,meas_noise_cov,meas_noise_mu,L)
% X_f = samples
% Y = observations
% Augmented Y = Ya
    yadd  = sum(X_f(1:2,:).*X_f(1:2,:),1);
    YP = [X_f(1:2,:)',yadd']; % Predicted measurements
    Ya = [Ya,L^2*ones(size(Ya,1),1)];
    meas_noise_mu =  [meas_noise_mu,0];
    meas_noise_cov = diag([meas_noise_cov,0.001]);
    W =  mvnpdf(Ya-YP,meas_noise_mu,meas_noise_cov);
end
