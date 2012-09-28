% DEMTOYPPCA1 Probabistic PCA demo on simulated data.
%
% FORMAT
% DESC
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011
%
% RCA

addpath(genpath('~/mlprojects/matlab/general/'))
importTool({'rca','ndlutil'})

% First we generate some artificial data.
% Observed, latent dimensionalities.
d = 3;
p = 2;
n = 4000;
% Low rank structure.
W = randn(d,p);
WWt = W*W';
% Noise variance.
sigma_n = .01;
% Generate observations from a low-rank covariance.
Theta = WWt + sigma_n*eye(d);
Y = gaussSamp(Theta, n);

% Plot data superimposed with Gaussian contours at 1 and 2 stds.
clf; hold on;
t = -pi:.02:pi;
k = length(t);
xyz = [zeros(k,1) cos(t)' sin(t)'];
[U,D] = eig(Theta);
A = (U*sqrt(D))'; % Level-set at 1 std.
for i = [1 2] 
    xyz_t = xyz*i*A;
    C = mvnpdf(xyz_t,0,Theta);
    scatter3(xyz_t(:,1), xyz_t(:,2), xyz_t(:,3), ones(k,1), C, '.');
end
C = mvnpdf(Y,0,Theta);
scatter3(Y(:,1), Y(:,2), Y(:,3), 20*ones(n,1), C,'.');
colormap('hot'), grid on; colorbar
