% DEMTOYEMRCA1 EM-RCA demo on simulated data.
%
% FORMAT
% DESC
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011
%
% RCA

clc
addpath(genpath('~/mlprojects/matlab/general/'))
addpath(genpath('~/mlprojects/rca/matlab/glasso/'))
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
importTool({'rca','ndlutil'})
asym = @(x) sum(sum(abs(x - x.'))); % Asymmetry test.

figure(1), clf, colormap('hot');    figure(2), clf, colormap('hot')
figure(3), clf, colormap('hot');    figure(4), clf, colormap('hot')
figure(5), clf, colormap('hot')

limit = 1e-4;
lambda = 5.^linspace(-8,3,30);
Sigma_hat = cell(length(lambda),1); Lambda_hat = cell(length(lambda),1);
triuLambda_hat = cell(length(lambda),1);
B = cell(length(lambda),1);
rocstats = zeros(length(lambda), 4);
options = struct('verbose',1,'order',-1);
emrca_options = struct('showProgress',0 , 'verbose',0, 'errorCheck',1);

%% Data generation.
d = 50; % Observed dimensions.
p = 3;  % Low-rank
n = 100;
sigma2_n = 1e-1; % Noise variance.
% Generate a sparse positive-definite precision matrix w/ given density and
% non-zero entries normally distributed w/ mean 1 and variance 2.
s = RandStream('mcg16807','Seed', 1985); RandStream.setDefaultStream(s) % 1985
density = .01;
spN = ceil((d^2-d)/2 * density);    % Number of non-zeros elements that satisfy the density.
sel = randperm(d^2);                % Random positions.
sel( mod(sel,d+1) == 1 ) = [];      % Remove positions of the diagonal.
sel = sel( 1:spN );                 % Enough random positions to satisfy the density.
Lambda = zeros(d);  Lambda(sel) = 1;        % Non-zero positions mask.
Lambda = triu(Lambda + Lambda',1);          % Upper triangular after balancing.
Lambda = Lambda .* (randn(d)*sqrt(2) + 1);  % Apply mask.
Lambda = Lambda + Lambda';          % Symmetrify.
pd = 1; diagonal = eye(d);
while pd > 0
    testLambda = Lambda + diagonal;
    [~, pd] = chol(testLambda);
    diagonal = diagonal*2;
end
Lambda = testLambda;    Sigma = pdinv(Lambda);
% Sample from p(Y).
W = randn(d,p);
WWt = W*W';
Theta = WWt + Sigma + sigma2_n*eye(d);
Y = gaussSamp(Theta, n);
Y = Y - repmat(mean(Y),n,1);
Cy = Y' * Y /n;
nonZero = find(ones(d));    % Induce any prior knowledge of zeros. If not, this is all ones.

figure(1), clf, colormap('hot')
subplot(131), imagesc(Lambda), title('sparse \Lambda'), colorbar
subplot(132), imagesc(Sigma), title('sparse-inverse \Sigma'), colorbar
subplot(133), imagesc(WWt), title('Low-rank WW'''), colorbar

%% Standard Glasso on (un)confounded simulated data, with varying lambda.
%{
confounders{1} = WWt;   confounders{2} = zeros(size(WWt));
for c = 1:2
    s = RandStream('mcg16807','Seed', 666); RandStream.setDefaultStream(s) % 23, 1e5
    Y_ = gaussSamp(confounders{c} + Sigma + sigma2_n*eye(d), n);   Y_ = Y_ - repmat(mean(Y_),n,1);
    Cy_ = Y_' * Y_/n;
    warmLambda_hat = pdinv(Cy_); % eye(d);
    funObj = @(x)sparsePrecisionObj(x, d, nonZero, Cy_);
    parfor (i = 1:length(lambda),8)
        Lambda_hat{i} = eye(d);
        Lambda_hat{i}(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda(i)*ones(d*d,1), options);
        figure(3), imagesc(Lambda_hat{i}), colormap(hot), colorbar, title([ 'GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
        rocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat{i});
    end
    TPs = rocstats(:,1); FPs = rocstats(:,2); FNs = rocstats(:,3); TNs = rocstats(:,4); 
    Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);    FPRs = FPs ./ (FPs + TNs);
    linestyle = {'-xb','--om'}; legends = {'GL','"Ideal" GL'};
    figure(2), hold on, plot(Recalls, Precisions, linestyle{c}), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
    figure(4), hold on, plot(FPRs, Recalls, linestyle{c}), xlim([0 1]), xlabel('FPR'), ylabel('TPR')
    AUCs(c) = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs); %#ok<SAGROW>
end
figure(2), legend(legends,1), title('Recall-Precision');
figure(4), legend([ legends{1} ' auc: ' num2str(AUCs(1)) ], [ legends{2} ' auc: ' num2str(AUCs(2)) ], 4), title('ROC');
%}

%% Recovery of low-rank component WWt by explaining away the true
% sparse-inverse covariance Sigma.
%{
Theta_exp = Sigma + sigma2_n*eye(d);
[S D] = eig(Cy, Theta_exp);    [D perm] = sort(diag(D),'descend');
W_hat = Theta_exp * S(:,perm(D>1)) * sqrt(diag(D(D>1)-1));
WWt_hat = W_hat * W_hat';
figure(3), clf, imagesc(WWt_hat), colorbar, title('WW'' by RCA');
figure(4), clf, imagesc(WWt_hat - WWt), title('WW''-WWt_hat'), colorbar;
%}

%% Recovery of sparse-inverse and low-rank covariance via iterative
% application of GLASSO and RCA.
sigma2_n = .01*trace(Cy); % 0.045 * trace(Cy)/d;        % Noise variance. (.045, .025)
[S D] = eig(Cy);     [D perm] = sort(diag(D),'descend');    % Initialise W with a PCA low-rank estimate.
W_hat_old = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));
WWt_hat_old = W_hat_old * W_hat_old';
Lambda_hat_old = pdinv(Cy); % eye(d);
tic
parfor (i = 1:length(lambda),8)    % Try different magnitudes of lambda.
    [WWt_hat_new, Lambda_hat_new, Lambda_hat_new_inv] = emrca(Y, WWt_hat_old, Lambda_hat_old, sigma2_n, lambda(i), nonZero, limit, emrca_options);
    % Plot results.
    figure(5), clf, colormap('hot')
    subplot(131), imagesc(Lambda_hat_new), colorbar, title([ 'GLasso/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    subplot(132), imagesc(Lambda_hat_new_inv), colorbar, title('\Sigma_{hat}'), colorbar
    WWt_hat_new(WWt_hat_new > max(WWt(:))) = max(WWt(:));   WWt_hat_new(WWt_hat_new < min(WWt(:))) = min(WWt(:));
    subplot(133), imagesc(WWt_hat_new), colorbar, title('RCA-recovered WW'''), colorbar
    % Performance stats. Row format in pstats : [ TP FP FN TN ].
    rocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat_new);
end
toc

%% Process performances measures.
TPs = rocstats(:,1); FPs = rocstats(:,2); FNs = rocstats(:,3); TNs = rocstats(:,4); 
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);

%% Plot performance.
figure(3), hold on, plot(Recalls, Precisions, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), title('Recall-Precision')
text(Recalls, Precisions, num2cell(lambda))
GLxy  = [1,0; 0.933,0; 0.866,0.025; 0.8,0.03; 0.8,0.04; 0.75,0.05;  0.6,0.06; 0.525,0.08; 0.4,0.1; 0.25,0.08; 0.25,0.1; 0.125,0.09; 0.05,0.125; 0.05,0.066; 0.05,1];
KGLxy = [1,0; 0.933,0; 0.8,0.01; 0.8,0.03; 0.733,0.04; 0.46,0.125; 0.41, 0.433; 0.2,0.6; 0.125,1];
IGLxy = [1,0; 0.933,0; 0.866,0.133; 0.6,0.2; 0.45,0.325; 0.125,0.666; 0.05,1];
plot(GLxy(:,1), GLxy(:,2), 'bo-', KGLxy(:,1), KGLxy(:,2), 'gx-', IGLxy(:,1), IGLxy(:,2), 'm.-')
legend('EM-RCA','Glasso','Kr-Glasso','Ideal Glasso')
figure(4), hold on, plot(FPRs, Recalls, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'RCA-GLasso auc: ' num2str(AUC) ], 4), title('ROC');
