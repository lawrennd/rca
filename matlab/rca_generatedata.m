function [Y,X,Z] = rca_generatedata(W,V,N,Mu,sigma_sq)

%   (Z)       (X)       | Y|Z,X ~ N(ΧW+ZV+μ, σ²I)
%     \       /         |   W,V ~ N(0,I)
%   V  \     /  W       | explained + residual covariance =
%       \   /           |----------------------------------
%        v v            |  ZZ' + σ² + XX' =
%  μ --> (Y) <-- σ²     |         Σ + ΧΧ'
% Generates data from a Residual Component Analysis model.
% Y = rca_generatedata(W,V,N,Mu,sigma_sq) generates N row vectors in Y
% distributed as illustrated in the graphical model above. The N rows of X
% and Z contain the N generated latent datapoints.
% author: Alfredo Kalaitzis, 17-11-09

% dimensionalities
[Q,D] = size(W);
P = size(V,1);
% true observation noise covariance
obs_noise = sigma_sq*eye(D);
% draw latent variables X,Z from zero-mean, unit-variance Gaussians.
% see 'help mgd_sample' for details
X = mgd_sample(zeros(N,Q), eye(Q));
Z = mgd_sample(zeros(N,P), eye(P));
% generate data from the true model
XW = X*W; ZV = Z*V;
Y = mgd_sample(XW + ZV + ones(N,1)*Mu, obs_noise);
end



