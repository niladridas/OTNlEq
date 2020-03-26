clc;clear;close;
X = [-20,-10,0,10,20; -8,-1,0,1,8];
% Data from the manifold : y = (x/10)^3
% We want the data to be in 1D manifold.

% The M matrix is constructed using Floyd's algorithm.
% There are total 5 data points
N = 5;
% Size of matrix M is 5x5
M_floyd  = [0, 12.21, 21.54, 31.59, 43.8; 
            12.21, 0, 10.05, 20.1, 32.31;
            22.26, 10.05, 0 , 10.05, 22.26;
            32.31, 20.1, 10.05, 0, 12.21;
            43.8, 31.59, 21.54, 12.21, 0];
% Constructing \tau(M)
S = M_floyd.^2;
H = eye(5,5)-(1/N)*ones(N,N);
tau_M = -H*S*H/2;

[V,D] = eig(tau_M);
D     = diag(D);

% We are interested in 1 dimensional representation
x_new = sqrt(D(1)).*V(:,1);
