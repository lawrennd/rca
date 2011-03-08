% specify mappings here.
W = [5 0 0;     % true W. Q-by-D.
    0 5 0];
V = [-3 4 5;    % true V. P-by-D.
    -1 -2 1];
N = 300;        % no. generated datapoints
Mu = [0 0 0];   % true bias
sigma_sq = .5;  % true observation noise variance

% observed and latent datapoints generated from RCA model
[Y,X,Z] = rca_generatedata(W,V,N,Mu,sigma_sq);
XW = X*W; ZV = Z*V;
% illustrate latent bases on Y space
figure(1), clf, grid on, axis equal, view([50 25]), hold on,
plot3([0 W(1,1)]*5, [0 W(1,2)]*5, [0 W(1,3)]*5, 'b->',...
    [0 W(2,1)]*5, [0 W(2,2)]*5, [0 W(2,3)]*5, 'b->',...
    [0 V(1,1)]*5, [0 V(1,2)]*5, [0 V(1,3)]*5, 'r->',...
    [0 V(2,1)]*5, [0 V(2,2)]*5, [0 V(2,3)]*5, 'r->','linewidth', 2),
% illustrate latent points on Y space
plot3(XW(:,1),XW(:,2),XW(:,3),'ob',...
    ZV(:,1),ZV(:,2),ZV(:,3),'or','markersize',4),
% illustrate observed data points
plot3(Y(:,1), Y(:,2), Y(:,3),'ok', 'markerfacecolor','g','markersize', 5),
hold off;
title('Y space'), legend('W:,1', 'W:,2', 'V:,1', 'V:,2', 'X', 'Z', 'Y');

% infer X through RCA model
[Xinf,D] = rca(Y, Z, sigma_sq);
Xinf = Xinf(:,1:2); % we only need the first 2 eigenvectors

% find latent rotation and scaling of inferred X
Xinf = Xinf*sqrt(N);
[d, Xtrns, tr] = procrustes(X, Xinf);
disp(['dissimilarity(sum of squared errors) = ' num2str(d)]);
% illustrate inferred X
figure(2), clf, axis equal,
subplot(1,2,1), hold on,
plot(X(:,1), X(:,2), 'rx'), legend('true X'),
subplot(1,2,2), hold on,
plot(Xtrns(:,1), Xtrns(:,2), 'gx'),legend('procrusted inferred X'),
for i = 1:size(X,1)
    subplot(1,2,1), text(X(i,1), X(i,2), num2str(i));
    subplot(1,2,2), text(Xtrns(i,1), Xtrns(i,2), num2str(i));
end
hold off;

