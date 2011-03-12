function [X,D] = rca(Y, Z, sigma_sq)
%{
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
% generalised problem Σ¯¹YY'S=SD, the solution being X as the first
% Q eigenvectors of YY'.
%}
%
% SEEALSO : eig
%
% Author: Alfredo A. Kalaitzis, 2009, 2011

% RCA

Sigma = Z*Z' + sigma_sq*eye(size(Z,1)); % Explained covariance.
%{
% Solve for S through a regular eigenvalue problem of transformed Y.
[U, Lamda] = eig(Sigma); % Eigen-decompose Sigma.
Lamda_sqrt = diag(sqrt(diag(Lamda)));
Lamda_invsqrt = diag(1./diag(Lamda_sqrt));
Y_ = Lamda_invsqrt*U'*Y; % Transform Y st K_ = X_*X_' + I
[V,D] = eig(Y_*Y_');
S = U*Lamda_invsqrt*V; % S = Σ¯¹Τ = UΛ¯¹U'UΛ¹/²V = UΛ¯¹/²V 
%}

Yinprod = Y*Y'; % Inner product matrix.
% Solve for S via a generalised eigenvalue problem of YY' and Σ: YY'S=ΣSD.
% S contains the eigenvectors of Yinprod. X columns are the Q principal
% eigenvectors, up to scaling and rotation.
[S,D] = eig(Yinprod,Sigma);
[D,ind] = sort(diag(D),'descend'); % Sort eigenvalues and save permutation.
X = S(:,ind); % Permute eigenvectors.
end
