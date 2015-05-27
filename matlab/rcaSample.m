function [Y,X,Z] = rcaSample(W,V,n,sigma_sq)

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

[d,q] = size(W); % Observed and latent dimensions.
p = size(V,2);
% Sample latent variables X,Z from zero-mean, unit-variance Gaussians.
X = gaussSamp(eye(q), n);
Z = gaussSamp(eye(p), n);
% XWt = X*W'; ZVt = Z*V';
Y = X*W' + Z*V' + gaussSamp(sigma_sq*eye(d), n); % Generate data from the model.
end



