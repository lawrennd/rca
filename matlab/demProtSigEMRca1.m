% DEMPROTSIGEMRCA1 EM-RCA demo on reconstruction of a protein
% signalling network.
%
% FORMAT
% DESC
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011, 2012
%
% RCA

clc
addpath(genpath('~/mlprojects/matlab/general/'))
importTool({'rca','ndlutil'})
addpath(genpath('~/mlprojects/rca/matlab/glasso/'))
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
addpath(genpath('~/mlprojects/rca/matlab/PROTSIG_DATA'))
% s = RandStream('mcg16807','Seed', 1985); RandStream.setDefaultStream(s) % 1985
% figure(1), clf, colormap('hot'); figure(2), clf, colormap('hot'); figure(3), clf, colormap('hot'); figure(4), clf, colormap('hot'); figure(5), clf, colormap('hot')

% General settings.
lambda = 5.^linspace(-8,3,50);
GLASSOrocstats = zeros(length(lambda), 4);
EMRCArocstats = zeros(length(lambda), 4);
options = struct('verbose',1,'order',-1);

%% Import data.
[Y1, colnames] = xlsread('1. cd3cd28.xls');
Y2 = xlsread('2. cd3cd28icam2.xls');
Y3 = xlsread('3. cd3cd28+aktinhib.xls');
Y = [Y1;Y2;Y3];
selection = randperm(size(Y,1));          % *** comment when doing stability selection ***
Y = Y(selection(1:fix(size(Y,1)/10)),:);    % Select 10% of the data. Comment when ALL data are used.

[n,d] = size(Y);
nonZero = find(ones(d));                    % Induce any prior knowledge of zeros. If not, this is all ones.
% Ground truth Lambda
Lambda = full(sparse([1,1,1, 2,2,2, 3,3, 4, 7, 8,8, 9,9, 1:11], [9,8,2, 9,8,6, 5,4, 5, 8, 11,10, 11,10, 1:11], ones(1,25)));
Lambda = Lambda + triu(Lambda,1)';
nodeLabels = {'praf','pmek','plcg','PIP2','PIP3','p44/42','pakts473','PKA','PKC','P38','pjnk'};

iSub = randsample(n,ceil(n*.9));            % Sub-sample from data for stability selection.
% Y = Y(iSub,:);                              % *** un-comment when doing stability selection ***

Y = Y - repmat(mean(Y), size(Y,1), 1);
Y = Y ./ repmat(std(Y), size(Y,1), 1);      % Normalise features.
Cy = Y' * Y / size(Y,1);
jit = randn(length(lambda),1)*.001;
% figure(1), colormap('hot'), imagesc(Lambda), title('Ground truth network'), colorbar, daspect('manual')

[B, U, jitter] = pdinv(Cy);
rocstats = emrcaRocStats(Lambda, B);        % Measure of how good B is for initializing Lambda_hat.

%% Standard Glasso on protein-signalling data, with varying lambda.
%
Lambda_hat_glasso = cell(length(lambda),1);
funObj = @(x)sparsePrecisionObj(x, d, nonZero, Cy);
warmLambda_hat = eye(d);
parfor (i = 1:length(lambda),8)
    Lambda_hat_glasso{i} = eye(d);
    Lambda_hat_glasso{i}(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda(i)*ones(d*d,1), options);
%     figure(3), imagesc(Lambda_hat_glasso{i}), colormap(hot), colorbar, title([ 'GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    % Evaluation
    GLASSOrocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat_glasso{i});
end
%}
%
TPs = GLASSOrocstats(:,1); FPs = GLASSOrocstats(:,2); FNs = GLASSOrocstats(:,3); % TNs = GLASSOrocstats(:,4); 
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
% FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);

figure(2), clf, colormap hot, hold on,
plot(Recalls, Precisions, '--xb'), %text(Recalls+jit, Precisions+jit, num2cell(lambda)),
xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision')

plot([1,.87,.27],[.275,.26,.57], 'g-', [1,.67,.2], [.275,.3,.5], 'b-')     % literature performance
% plot(rocstats(1)./ (rocstats(1) + rocstats(3)), rocstats(1) / (rocstats(1) + rocstats(2)), 'rs') % Empirical inv.cov performance.
legend('Glasso','Kronecker-Glasso','Glasso (reported)', 'inv.cov')
% figure(4), clf, hold on, plot(FPRs, Recalls, '-xb'), xlim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'GL auc: ' num2str(AUC) ], 4), title('ROC');
%}

%% Recovery of sparse-inverse and low-rank covariance via iterative application of GLASSO and RCA.
lambda = 5.^linspace(-8,3,50);
emrca_options = struct('limit',1e-4, 'showProgress',0 , 'verbose',0, 'errorCheck',1, 'maxNumIter',1000);
sigma2_n =  0.5*trace(Cy)/d; % 0.3 * trace(Cy)/d;        % Noise variance. (0.3)
[S D] = eig(Cy);     [D perm] = sort(diag(D),'descend');    % Initialise W with a PCA low-rank estimate.
W_hat_old = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));
WWt_hat_old = W_hat_old * W_hat_old';
Lambda_hat_emrca = cell(length(lambda),1);
Lambda_hat_old = eye(d);
tic
parfor (i = 1:length(lambda),8)            % Try different magnitudes of lambda.
    [WWt_hat_new, Lambda_hat_emrca{i}] = emrca(Y, WWt_hat_old, Lambda_hat_old, sigma2_n, lambda(i), nonZero, emrca_options);
    % Plot results.
    %{
    figure(5), clf, subplot(131), imagesc(Lambda_hat_emrca{i}), colorbar, title([ 'EM/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    % subplot(132), imagesc(Lambda_hat_new_inv), colorbar, title('\Sigma_{hat}'), colorbar
    subplot(133), imagesc(WWt_hat_new), colorbar, title('RCA-recovered WW'''), colorbar
    %}
    EMRCArocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat_emrca{i});     % Performance stats. Row format in pstats : [ TP FP FN TN ].
end
toc

%% Process performances measures.
%{
TPs = EMRCArocstats(:,1); FPs = EMRCArocstats(:,2); FNs = EMRCArocstats(:,3); % TNs = EMRCArocstats(:,4); 
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);

%% Plot performance.
figure(2), hold on, plot(Recalls, Precisions, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
text(Recalls+jit, Precisions+jit, num2cell(lambda))
plot([1,.87,.27],[.275,.26,.57], 'g-', [1,.67,.2], [.275,.3,.5], 'b-')   % Literature performance.
legend('Glasso','EM-RCA','Kronecker-Glasso (reported)','Glasso (reported)')

% figure(4), hold on, plot(FPRs, Recalls, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'EM/RCA auc: ' num2str(AUC) ], 4), title('ROC');
%}