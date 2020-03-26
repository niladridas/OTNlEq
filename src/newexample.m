clc;close;clear;
format longE;
% Example taken from the paper: http://www-personal.umich.edu/~dsbaero/library/ConferencePapers/ACC2008/UKFEquConstraints.pdf
% ATITUDE ESTIMATION
Tsample = 0.1;
%% Simulate the actual attitudedyn system
e0      = [0.9603 0.1387 0.1981 0.1387];
u       = @(t) [0.03*sin(2*pi*t/600) 0.03*sin((2*pi/600)*t-300) 0.03*sin((2*pi/600)*t-600)];
Tend    = 1000;
% tspan   = 0:Tsample:Tend;
options = odeset('RelTol',1e-12,'AbsTol',1e-10);
% [T,Econt]   = ode45(@(t,e) attitudedyn(t,e,u),tspan,e0,options);
% norm2Econt = sqrt(sum(Econt.*Econt,2));
% semilogy((1-norm2Econt))

%% Simulate the discretized attitudedyn system
beta = [0.001 -0.001 0.0005]';
Qm = 10^(-5)*eye(3);
Qb = 10^(-10)*eye(3);
Nsamples = (Tend/Tsample);
e = e0';
Edisc = zeros(Nsamples+1,4);
Edisc(1,:) = e;
for k = 1:Nsamples
    [e,beta] = attitudediscretedyn(k,e,beta,Qb,Tsample,u);
    Edisc(k+1,:) = e';
end

%% Percent error in each of the quaternions between continuous and discrete dynamics
% pE = sqrt((1/size(Econt,1))*sum(((Econt-Edisc)./Econt).^2,1));

%% Measurement of attitude measuring system
Y_meas = zeros(size(Edisc,1),6);
R = 10^(-4)*eye(6);
for i = 1:size(Edisc,1)
    Y_meas(i,:) = discmeasmodel(Edisc(i,:),R);
end    
plot(Y_meas)
%% The measurement Model runs at at much slower rate of 10*Tsample
% For this we extract the measurements as
Y_lowratemeas = Y_meas(1:10:end,:);
Nmeas = size(Y_lowratemeas,1);
%% OT filtering without constraints
Ns = 20; % sample size
P0 = 0.05*eye(4); % Initial covariance which is different from the paper
% Generate initial samples satisfying the constraint.
X0_pert = mvnrnd(zeros(1,4),P0,Ns);
X0 = repmat(e0',1,20)+X0_pert';
%%
% OT filtering
X_prior = X0;
cost    = @(X) distance_matrix(X);
weight  = @(X,Y) weightPM(X,Y,R);
X_mean  = zeros(Nmeas,4);
X_cov   = zeros(Nmeas,4); % sig_x1x1, sig_x2x2, sig_x3x3, sig_x4x4

%%
for i = 1:Nmeas
    % UPDATE
    X_post = OT_filter(X_prior,Y_lowratemeas(i,:),cost,weight,@OT_constants,@Optimal_Transport);
%     X_post = OT_filtertnonlineq(X_prior,Y_meas(i,:),cost,weight,@OT_constants,@Optimal_Transport,L);
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
    if i==Nmeas
        break;
    end
    tspan    = [(i-1)*10*Tsample,i*10*Tsample];
    for j = 1:Ns
        [~,Xtmp]   = ode45(@(t,e) attitudedyn(t,e,u),tspan,X_post(:,j),options);
        X_prior(:,j) =  Xtmp(end,:)';
    end
end
%%
plot(sqrt(sum(X_mean.*X_mean,2)))

%%
function W = weightPM(E,Ya,R)
% X_f = samples
% Y = observations
    Ns = size(E,2);
    W = zeros(Ns,1);
    for i = 1:Ns
        b = E(1,i); c = E(2,i); d = E(3,i); a = E(4,i);
        C = [2*a^2-1+2*b^2, 2*(b*c+a*d), 2*(b*d-a*c);
             2*(b*c-a*d), 2*(a^2+c^2)-1, 2*(c*d+a*b);
             2*(b*d+a*c), 2*(c*d-a*b), 2*(a^2+d^2)-1];
        y1 = C*[1;0;0];
        y2 = C*[0;1;0];
        YP = [y1;y2];
        W(i,1) =  mvnpdf(Ya-YP',zeros(1,6),R);
    end
end

