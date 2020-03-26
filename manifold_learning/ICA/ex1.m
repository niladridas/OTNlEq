% Author: Niladri Das
% Date  : 11th Oct 2019
% Paper : Spectral independent component analysis - Amit Singer
clear;close;clc;
n = 10;                    % Denotes independent random variables.
A = RandOrthMat(n);        % Generate orthogonal matrix
N = 40;                    % Number of mixed data generated

rng(123);
mu = zeros(1,n); sig = ones(n,n);
S = mvnrnd(mu, sig, N);
S = S';                    % Actual source variable
X = A*S;                   % Observed variable

% Calculate the covariance of the data
C     = cov(X');
[V,D] = eig(C);
D     = diag(D);

W = gen_Wmatrix(X);
% Normalizing W by D. D is a diagonal matrix with diagonals having row sum
% of W.
D = diag(sum(W,2));
% The graph laplacian is:
L = diag(1./sum(W,2))*W - eye(size(W,1));
% The top few eigen values and eigen vector of this graph laplacian is used
% for dimensionality reduction or clustering
%% When X is samples uniformly and independently from manifold M, the graph laplacian 
%  converges to the Laplace-Beltrami operator of the manifold.
%%
% The joint differential entrophy H(S) of independent random variable = Sum 
% of their marginal entrophies Sum H(S_j).
% Joint differential entrophy is invariant to orthogonal transformations;
% H(AS) = H(S), if A is orthogonal
%
% Optimization problem: Find A such that Sum H(S_j) is max or Sum
% H((A^-1*X)_j) is max.
% Global search algorithm: Newton or a gradient method.
% Disadvantages: Get stuck in a local maxima.
% This technique based on differential entrophy.
%%
% There is another technique based on graph Laplacian of data.
% Graph Laplacian and the backward Focker Plank operator
function [W] = gen_Wmatrix(X)
    % X = dimension of X times number of samples
    N = size(X,2);
    W = zeros(N,N);
    % Calculate the bandwidth
    % Median based bandwidth is implemented here:
    % 1. Calculate the Euclidean distance for all pairs of points
    D = zeros(N,N);
    for i = 1:(N-1)
        for j = (i+1):N
            D(i,j) = dot(X(:,i),X(:,j));
            D(j,i) = D(i,j);
        end
    end
    epsilon = median(vec(D(D>0)));
    kernel = @(x,y) exp(-(norm(x-y).^2)/(2*epsilon));
    %
    for i = 1:(N-1)
        for j = (i+1):N
            W(i,j) = kernel(X(:,i),X(:,j));
            W(j,i) = W(i,j);
        end
    end
end