clc;clear;close;
X = [-20,-10,0,10,20; -8,-1,0,1,8];
% Data from the manifold : y = (x/10)^3
% We want the data to be in 1D manifold.

% Step 1:
% Calculate the covariance of the data
C = cov(X');
[V,D] = eig(C);

% Step 2:
% One eigen value is comparatively larger than the other
% The principle axis is:
paxis = V(:,2)/norm(V(:,2));

X_new = paxis'*X;

% Plot of all points:
x1 = linspace(-20,20,100);
plot(x1,(x1/10).^3); hold on;
scatter(X(1,:),X(2,:),'filled','k'); hold on;
plot(x1,(paxis(2)/paxis(1))*x1); hold on;
X_proj = (paxis(1)/sqrt(paxis(1)^2+paxis(2)^2))*X_new;
Y_proj = (paxis(2)/sqrt(paxis(1)^2+paxis(2)^2))*X_new;
scatter(X_proj,Y_proj,'filled','red'); hold on;
grid on;

