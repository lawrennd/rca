% DEMTOYEMRCA1 EM-RCA demo on simulated data.
%
% FORMAT
% DESC
%
% SEEALSO : emrca
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011, 2012
%
% RCA

% clc
addpath(genpath('~/mlprojects/matlab/general/'))
addpath(genpath('~/mlprojects/rca/matlab/glasso/'))
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
importTool({'rca','ndlutil'})
asym = @(x) sum(sum(abs(x - x.'))); % Asymmetry test.
% figure(1), clf, colormap('hot'); figure(2), clf, colormap('hot'); figure(3), clf, colormap('hot'); figure(4), clf, colormap('hot'); figure(5), clf, colormap('hot')

% General settings.
lambda = 5.^linspace(-8,3,30);
rocstats = zeros(length(lambda), 4);
GLASSOrocstats = {zeros(length(lambda), 4) zeros(length(lambda), 4)};
EMRCArocstats = zeros(length(lambda), 4);
options = struct('verbose',1,'order',-1);

%% Synthetic data.
d = 20; % Observed dimensions.
p = 3;  % Low-rank
n = 100;
density = .01;
%{
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
%}
load simdata_Y.mat;
%{
    %     s = RandStream('mcg16807','Seed', 1e+5); RandStream.setDefaultStream(s) % 1985
    [Lambda, Sigma] = genspd(d, density, 1, 2); % Generate a sparse positive-definite precision matrix.
    X = randn(n,p);
    W = randn(p,d);
    %     W = reshape(gaussSamp(kron(Sigma,eye(p)),1), p,d);
    XW = X*W;
    V = gaussSamp(Sigma, n);
    %     V = reshape(gaussSamp(kron(Sigma,eye(n)),1), n,d);
    E = randn(n,d);
 %}    
rho = sqrt((trace(XW'*XW)/trace(V'*V)) * 1);
sigma_n = sqrt(((trace(XW'*XW) + trace((rho^2)*(V'*V)))/trace(E'*E)) * .1);
Y = XW + rho*V + sigma_n*E;
Y_unc = rho*V + sigma_n*E;
%{
% Sample from p(Y).
W = randn(d,p);
WWt = W*W';
Theta = WWt + Sigma + sigma2_n*eye(d);
Y = gaussSamp(Theta, n);
%}

iSub = randsample(n,ceil(n*.9));        % Sub-sample from data.
Y = Y(iSub,:);                          % *** un-comment when doing stability selection ***
Y_unc = Y_unc(iSub,:);

Y = Y - repmat(mean(Y),size(Y,1),1);
Y_unc = Y_unc - repmat(mean(Y_unc),size(Y_unc,1),1);
Cy = Y' * Y / size(Y,1);
nonZero = find(ones(d));                % Induce any prior knowledge of zeros. If not, this is all ones.

%{
figure(1), subplot(131), imagesc(Lambda), title('sparse \Lambda'), colorbar
subplot(132), imagesc(Sigma), title('sparse-inverse \Sigma'), colorbar
subplot(133), imagesc(X*X'), title('Low-rank XX'''), colorbar
%}

%% Standard Glasso on (un)confounded simulated data, with varying lambda.
%{
Y_ = {Y, Y_unc};
linestyle = {'-xb','--om'}; legends = {'GL','"Ideal" GL'};
Lambda_hat_glasso = cell(length(lambda),2);
for c = 1:2
    Cy_ = Y_{c}'*Y_{c} / size(Y_{c},1);
    warmLambda_hat = eye(d);
    funObj = @(x)sparsePrecisionObj(x, d, nonZero, Cy_);
    parfor (i = 1:length(lambda),8)
        Lambda_hat_glasso{i,c} = eye(d);
        Lambda_hat_glasso{i,c}(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda(i)*ones(d*d,1), options);
        figure(3), imagesc(Lambda_hat_glasso{i,c}), colormap(hot), colorbar, title([ 'GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
        rocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat_glasso{i,c});
    end
    GLASSOrocstats{c} = rocstats;
    TPs = GLASSOrocstats{c}(:,1); FPs = GLASSOrocstats{c}(:,2); FNs = GLASSOrocstats{c}(:,3); % TNs = GLASSOrocstats{c}(:,4); 
    Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);   % FPRs = FPs ./ (FPs + TNs);
    figure(2), hold on, plot(Recalls, Precisions, linestyle{c}), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
    %     figure(4), hold on, plot(FPRs, Recalls, linestyle{c}), xlim([0 1]), xlabel('FPR'), ylabel('TPR')
    %     AUCs(c) = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);
end
figure(2), legend(legends,1)
% figure(4), legend([ legends{1} ' auc: ' num2str(AUCs(1)) ], [ legends{2} ' auc: ' num2str(AUCs(2)) ], 4), title('ROC');
%}

%% Recovery of sparse-inverse and low-rank covariance via iterative application of EM and RCA.
emrca_options = struct('limit',1e-4, 'showProgress',0 , 'verbose',0, 'errorCheck',1, 'maxNumIter',1000);
sigma2_n = 0.5*trace(Cy)/d;          % 0.045 * trace(Cy)/d;        % Noise variance. (.045)
[S D] = eig(Cy);     [D perm] = sort(diag(D),'descend');    % Initialise W with a PCA low-rank estimate.
W_hat_old = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));
WWt_hat_old = W_hat_old * W_hat_old';
Lambda_hat_emrca = cell(length(lambda),1);
Lambda_hat_old = eye(d);
tic
parfor (i = 1:length(lambda),8)     % Try different magnitudes of lambda.
    [WWt_hat_new, Lambda_hat_emrca{i}, Lambda_hat_new_inv] = emrca(Y, WWt_hat_old, Lambda_hat_old, sigma2_n, lambda(i), nonZero, emrca_options);
    % Plot results.
    %{
    figure(5), clf, colormap('hot')
    subplot(131), imagesc(Lambda_hat_emrca{i}), colorbar, title([ 'EM/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    subplot(132), imagesc(Lambda_hat_new_inv), colorbar, title('\Sigma_{hat}'), colorbar
%     WWt_hat_new(WWt_hat_new > max(WWt(:))) = max(WWt(:));   WWt_hat_new(WWt_hat_new < min(WWt(:))) = min(WWt(:));
%     subplot(133), imagesc(WWt_hat_new), colorbar, title('RCA-recovered WW'''), colorbar
    %}
    EMRCArocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat_emrca{i});    % Performance stats. Row format in pstats : [ TP FP FN TN ].
end
toc

%% Process performances measures.
%{
TPs = EMRCArocstats(:,1); FPs = EMRCArocstats(:,2); FNs = EMRCArocstats(:,3); TNs = EMRCArocstats(:,4); 
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);

%% Plot performance.
figure(2), hold on, plot(Recalls, Precisions, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
% text(Recalls, Precisions, num2cell(lambda))
% GLxy  = [1,0; 0.933,0; 0.866,0.025; 0.8,0.03; 0.8,0.04; 0.75,0.05;  0.6,0.06; 0.525,0.08; 0.4,0.1; 0.25,0.08; 0.25,0.1; 0.125,0.09; 0.05,0.125; 0.05,0.066; 0.05,1];
% KGLxy = [1,0; 0.933,0; 0.8,0.01; 0.8,0.03; 0.733,0.04; 0.46,0.125; 0.41, 0.433; 0.2,0.6; 0.125,1];
% IGLxy = [1,0; 0.933,0; 0.866,0.133; 0.6,0.2; 0.45,0.325; 0.125,0.666; 0.05,1];
% plot(GLxy(:,1), GLxy(:,2), 'bo-', KGLxy(:,1), KGLxy(:,2), 'gx-', IGLxy(:,1), IGLxy(:,2), 'm.-')
legend('GL','Ideal GL','EM-RCA')
% figure(4), hold on, plot(FPRs, Recalls, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'EM/RCA auc: ' num2str(AUC) ], 4), title('ROC');
%}