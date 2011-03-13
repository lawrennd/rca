function [X,D] = rca(Y, varargin)
%
%   (Z)       (X)       | Y|Z,X ~ N(ΧW'+ZV'+μ, σ²I)
%     \       /         |   W,V ~ N(0,I)
%   V  \     /  W       |-----------------------------------
%       \   /           | residual + explained covariance
%        v v            |    XX'   +      ZZ' + σ² 
%  μ --> (Y) <-- σ²     |    ΧΧ'   +          Σ
%
% Problem: We consider the case where Z is known and hence the covariance
% of data Y is already partially explained by Σ=ZZ'+σ². One would want to
% estimate X by marginalising out the mapping W and optimise the observed
% data marginal N(Y|0,XΣwX'+σ²I) wrt X. The problem is to separate the X
% from Σw after estimation of the product XΣwX'.
%
% Spectral solution: If we transform the observations Y st X'X = I
% (equivalent to constraining the prior over X to be spherical) we
% eliminate the aforementioned indeterminancy between X and Σw while at the
% same time we maintain the global optima in the solution of the
% generalised problem YY'S=ΣSD, the solution being X as the first
% Q generalised eigenvectors of YY' and Σ.
%
% SEEALSO : eig
%
% Author: Alfredo A. Kalaitzis, 2009, 2011

% RCA

if length(varargin) > 1
    Z = varargin{1};
    sigma_sq = varargin{2};
    Sigma = Z*Z' + sigma_sq*eye(size(Z,1)); % Explained covariance.
else
    Sigma = varargin{1};
end
%{
% Solve for S through a regular eigenvalue problem of transformed Y.
[U, Lamda] = eig(Sigma); % Eigen-decompose Sigma.
Lamda_sqrt = diag(sqrt(diag(Lamda)));
Lamda_invsqrt = diag(1./diag(Lamda_sqrt));
Y_ = Lamda_invsqrt*U'*Y; % Transform Y st K_ = X_*X_' + I
[V,D] = eig(Y_*Y_');
S = U*Lamda_invsqrt*V; % S = Σ¯¹Τ = UΛ¯¹U'UΛ¹/²V = UΛ¯¹/²V 
%}

YYt = Y*Y'; % Inner product matrix; covariance in the dual-space of Y.
% Solve for S via a generalised eigenvalue problem of YY' and Σ: YY'S=ΣSD.
% S contains the generarised eigenvectors of YY' and Σ. X columns are the Q
% principal generalised eigenvectors, up to scaling and rotation.
[S,D] = eig(YYt,Sigma);
[D,ind] = sort(diag(D),'descend'); % Sort eigenvalues and save permutation.
X = S(:,ind); % Permute eigenvectors.
% T = Sigma*S;
end
