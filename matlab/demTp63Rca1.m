% DEMTP63RCA1 RCA demo on TP63 expression time-series.
% FORMAT
% DESC Fit residual components (row rank covariance), with RCA, which
% potentially explain differences of gene expression between case and
% control. The explained covariance is an RBF kernel which encodes for
% the two conditions having the same behaviour. Tests on Della Gatta data.
%
% COPYRIGHT: Alfredo A. Kalaitzis, 2011
%
% SEEALSO : kernCreate, kernCompute, compareROC
%
% RCA

clear
addpath(genpath('~/mlprojects/matlab/general/'))
importLatest('netlab')
importTool({'kern','optimi','datasets','gprege'})

load DellaGattaData.mat
N = size(exprs_tp63_RMA,2);

tTrue = [timepoints; truncTimepoints]; % [0:20:240 0 20 40 60 120 180 240]';
kern = kernCreate(tTrue, 'rbf');
typDiff = 20; % Typical time difference between samples.
kern.inverseWidth = 1/(typDiff*typDiff); % Characteristic inverse-width.

% First kernel assumes the two profiles are generated by the same process;
% The second, that the two profiles are generated independently.
noise = 1e-4;
K_comb = kernCompute(kern, tTrue) + eye(20)*noise;

    % for k=0:20;
    % idx_top = sortIndex( ((k*1000)+1) : ((k*1000)+1000) );
    % idx_top = sortIndex(1:500);
idx_top = 1:length(exprs_tp63_RMA);

% Concatenate profiles of both conditions.
Y = [exprs_tp63_RMA(:, idx_top); exprs_null_RMA(:, idx_top)];

% In this application we apply the dual representation of RCA. This mean
% that want centred data in the sense that the mean is zero across the
% features, not the samples (take off the mean of what is larger). This
% means that in this case the mean is computed across the genes, not the
% timepoints. Here Y is sized Samples x Genes.
Y = Y - repmat(mean(Y,2),1, size(Y,2));

% We usually normalise by the number of samples, but in the dual case, we
% normalise by the number of features.
Cy = Y*Y'/size(Y,2);

[S D] = eig(Cy, K_comb);    [D perm] = sort(diag(D),'descend');
X = K_comb * S(:,perm(D>1)) * sqrt(diag(D(D>1)-1)); % Retrieve residual variance basis.
    % [X D] = rca(Y', K_comb);

% Plot evectors, combined kernel.
figure(1), clf, colormap gray
plot(X(:,1:end)), xlim([0 21]), title('Generalised eigenvectors (combined case)')
imagesc(K_comb), colorbar, daspect([1 1 1])

figure(2), clf
Y = [exprs_tp63_RMA; exprs_null_RMA]; % Prepare the whole dataset.
Y = Y - repmat(mean(Y,2),1, N); % Remove mean across features (genes).
Xproj = S(:,perm(D>1))'*Y;
dists = sum(Xproj.^2,1); % Compute norms.

% auc = compareROC(dists, DGatta_labels_byTSNItop100);
auc = compareROC(dists, DGatta_labels_byTSNItop100, BATSrankingFilepath);

    