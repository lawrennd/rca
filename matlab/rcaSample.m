function [Y,X,Z] = rcaSample(W,V,N,sigma_sq)

%   (Z)       (X)       | Y|Z,X ~ N(ΧW'+ZV'+μ, σ²I)
%     \       /         |   W,V ~ N(0,I)
%   V  \     /  W       |-----------------------------------
%       \   /           | residual + explained covariance
%        v v            |=   XX'   + ZZ' + σ² 
%  μ --> (Y) <-- σ²     |=   ΧΧ'   +     Σ
%
% Generates data from a RCA (residual components analysis) model.
% Y = rcaSample(W,V,N,sigma_sq) generates N row vectors in Y
% distributed as illustrated in the graphical model above. The N rows of X
% and Z contain the N generated latent datapoints.
%
% SEEALSO : gaussSamp
%
% Author: Alfredo A. Kalaitzis, 2009, 2011

% RCA

[D,Q] = size(W); % Observed and latent dimensions.
P = size(V,2);
% Sample latent variables X,Z from zero-mean, unit-variance Gaussians.
X = gaussSamp(eye(Q), N);
Z = gaussSamp(eye(P), N);
% XWt = X*W'; ZVt = Z*V';
Y = X*W' + Z*V' + gaussSamp(sigma_sq*eye(D), N); % Generate data from the model.
end



