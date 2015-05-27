% DEMRCA Test RCA on low-dim toy data.
%
% SEEALSO : rcaSample, rca, procrustes
%
% Author: Alfredo A. Kalaitzis, 2009, 2011

% RCA

addpath(genpath('~/mlprojects/matlab/'))
addpath(genpath('~/mlprojects/gp/matlab/'));
addpath(genpath('~/mlprojects/ndlutil'));
addpath(genpath('~/mlprojects/rca/matlab/'));


% Specify mappings here.
W = [5 0 0; 0 5 0]';    % True W. D-by-Q.
V = [-3 4 5; -1 -2 1]'; % True V. D-by-P.
N = 100;                % No. generated datapoints.

sigma_sq = .01;  % True noise variance.

% figure(4), clf
c = 1;

for sigma_sq = logspace(-2, 2, 20)
    
randn('seed', 1e7)
rand('seed', 1e7)

% Sample latent data-points and generate the observed.
[Y, X, Z] = rcaSample(W, V, N, sigma_sq);
XWt = X*W'; ZVt = Z*V';

% Illustrate latent bases in Y space.
z = zeros(2,1); W=W*4; V=V*4;
figure(3), clf, grid on, axis equal, hold on, view([50 25])
xlabel('X'), ylabel('Y'), zlabel('Z')
quiver3(z, z, z, W(1,:)', W(2,:)', W(3,:)', 1, 'b', 'linewidth', 2);
quiver3(z, z, z, V(1,:)', V(2,:)', V(3,:)', 1, 'r', 'linewidth', 2);

% Illustrate latent points on Y space.
plot3(  XWt(:,1),XWt(:,2),XWt(:,3), 'ob', 'markersize', 5, 'linewidth',1);
plot3(  ZVt(:,1),ZVt(:,2),ZVt(:,3), 'or', 'markerfacecolor','r','markersize', 5);
% Illustrate observed data points.
plot3(Y(:,1), Y(:,2), Y(:,3),'ok', 'markerfacecolor','g','markersize', 6);
title('Y space');


fill3(  [W(1,1)+W(1,2) W(1,1)-W(1,2) -W(1,1)-W(1,2) -W(1,1)+W(1,2)], ...
        [W(2,1)+W(2,2) W(2,1)-W(2,2) -W(2,1)-W(2,2) -W(2,1)+W(2,2)], ...
        [W(3,1)+W(3,2) W(3,1)-W(3,2) -W(3,1)-W(3,2) -W(3,1)+W(3,2)], ...
        'b', 'FaceAlpha', .2, 'EdgeColor', [0 0 1], 'EdgeAlpha', 0.1);
fill3(  [V(1,1)+V(1,2) V(1,1)-V(1,2) -V(1,1)-V(1,2) -V(1,1)+V(1,2)], ...
        [V(2,1)+V(2,2) V(2,1)-V(2,2) -V(2,1)-V(2,2) -V(2,1)+V(2,2)], ...
        [V(3,1)+V(3,2) V(3,1)-V(3,2) -V(3,1)-V(3,2) -V(3,1)+V(3,2)], ...
        'r', 'FaceAlpha', .2, 'EdgeColor', [1 0 0], 'EdgeAlpha', 0.1);
set(gcf,'renderer','painters')


%% Recover X through RCA.
[Xinf,D] = rca(Y', Z, sigma_sq);
Xinf = Xinf(:,1:2); % Need only the first 2 eigenvectors.
% Find rotation and scaling of recovered X.
Xinf = Xinf*sqrt(N);
Xtrns = Xinf;
[~, Xtrns, tr] = procrustes(X, Xinf);
d(c) = sum(sum((Xtrns-X).^2,2))/N; %sum(sum((X-repmat(mean(X,1),size(X,1),1)).^2,1));
c = c+1;

% disp(['Dissimilarity (sum of squared errors) = ' num2str(d)]);

%% Plot inferred X.
numstrs = num2cell((1:size(X,1)));
figure(2), clf,
% subplot(1,2,1),
hold on, plot(X(:,1), X(:,2), 'xg', 'MarkerSize', 3),
text(X(:,1), X(:,2), numstrs, 'color', [0 .7 0]),
hold on, plot(Xtrns(:,1), Xtrns(:,2), 'xr', 'MarkerSize', 3)
text(Xtrns(:,1), Xtrns(:,2), numstrs, 'color', [.8 0 0])
legend('True','Recovered')
% title('True')
% subplot(1,2,2),
% title('Inferred (up to rotation/scaling)')
end

figure(4), hold on, loglog(logspace(-2, 2, 20), d,'x-', 'color', [.8 0 0]), xlabel('\sigma^2'), ylabel('MSE')

