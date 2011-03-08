function [X,D] = rca(Y, Z, sigma_sq)

%   (Z)       (X)       | Y|Z,X ~ N(ΧW+ZV+μ, σ²I)
%     \       /         |   W,V ~ N(0,I)
%   V  \     /  W       | explained + residual covariance =
%       \   /           |----------------------------------
%        v v            |  ZZ' + σ² + XX' =
%  μ --> (Y) <-- σ²     |         Σ + ΧΧ'
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
% Author: Alfredo Kalaitzis, 18-11-09

Sigma = Z*Z' + sigma_sq*eye(size(Z,1)); % explained covariance
% Solve for S through a regular eigenvalue problem of transformed Y.
% [U, Lamda_sq] = eig(Sigma);
% Lamda = diag(sqrt(diag(Lamda_sq)));
% Lamda_inv = diag(1./diag(Lamda));
% Yprm = Lamda_inv*U'*Y; % transform Y st Kprm = Xprm*Xprm' + I
% [V,D] = eig(Yprm*Yprm');
% S = U*Lamda_inv*V; % S = Σ¯¹Τ = UΛ¯²U'UΛV = UΛ¯¹V

Yinprod = Y*Y'; % inner product matrix
% solve for S via a generalised eigenvalue problem of Y.
[S,D] = eig(Yinprod,Sigma);

% S now contains the eigenvectors of Yinprod and thus X as the Q principal
% eigenvectors.
[D,ind] = sort(diag(D),'descend'); % sort eigenvalues and keep permutation
X = S(:,ind); % sort eigenvectors according to permutation
end
