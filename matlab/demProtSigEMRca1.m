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
addpath(genpath('~/mlprojects/rca/matlab/glasso/'))
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
addpath(genpath('~/mlprojects/rca/matlab/PROTSIG_DATA'))
importTool({'rca','ndlutil'})
s = RandStream('mcg16807','Seed', 1985); RandStream.setDefaultStream(s) % 1985

figure(1), clf, colormap('hot');    figure(2), clf, colormap('hot')
figure(3), clf, colormap('hot');    figure(4), clf, colormap('hot')
figure(5), clf, colormap('hot')

limit = 1e-4;
lambda = 5.^linspace(-8,3,50);
Sigma_hat = cell(length(lambda),1); Lambda_hat = cell(length(lambda),1);
triuLambda_hat = cell(length(lambda),1);
B = cell(length(lambda),1);
rocstats = zeros(length(lambda), 4);
emrca_options = struct('showProgress',0 , 'verbose',0, 'errorCheck',1);
options = struct('verbose',1,'order',-1);


%% Import data.
[Y1, colnames] = xlsread('1. cd3cd28.xls');
Y2 = xlsread('2. cd3cd28icam2.xls');
Y3 = xlsread('3. cd3cd28+aktinhib.xls');
Y = [Y1;Y2;Y3];
selection = randperm(size(Y,1));
Y = Y(selection(1:fix(size(Y,1)/10)),:);    % Select 10% of the data.
[n,d] = size(Y);
nonZero = find(ones(d));                    % Induce any prior knowledge of zeros. If not, this is all ones.
% Ground truth Lambda
Lambda = full(sparse([1,1,1, 2,2,2, 3,3, 4, 7, 8,8, 9,9, 1:11], [9,8,2, 9,8,6, 5,4, 5, 8, 11,10, 11,10, 1:11], ones(1,25)));
Lambda = Lambda + triu(Lambda,1)';
Y = Y - repmat(mean(Y),n,1);
Y = Y ./ repmat(std(Y),n,1);                % Normalise features.
Cy = Y' * Y / n;
funObj = @(x)sparsePrecisionObj(x, d, nonZero, Cy);
warmLambda_hat = pdinv(Cy); %eye(d);
jit = randn(length(lambda),1)*.001;
figure(1), colormap('hot'), imagesc(Lambda), title('Ground truth network'), colorbar, daspect('manual')

%% Standard Glasso on protein-signalling data, with varying lambda.
%{
A = boolean( triu(Lambda,1) ~= 0 );
for i = 1:length(lambda)
    Lambda_hat{i} = eye(d);
    Lambda_hat{i}(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda(i)*ones(d*d,1), options);
    triuLambda_hat{i} = triu(Lambda_hat{i}, 1);
    figure(3), imagesc(Lambda_hat{i}), colormap(hot), colorbar, title([ 'GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    % Evaluation
    rocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat{i});
end
TPs = rocstats(:,1); FPs = rocstats(:,2); FNs = rocstats(:,3); TNs = rocstats(:,4); 
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);

figure(2), clf, hold on, plot(Recalls, Precisions, '-xb'), text(Recalls+jit, Precisions+jit, num2cell(lambda)),
xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), plot([1,.87,.27],[.275,.26,.57], 'g-', [1,.67,.2], [.275,.3,.5], 'b-')     % literature performance
legend('Glasso','Kronecker-Glasso','Glasso (reported)'), title('Recall-Precision');
figure(4), clf, hold on, plot(FPRs, Recalls, '-xb'), xlim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'GL auc: ' num2str(AUC) ], 4), title('ROC');
%}

%% Recovery of sparse-inverse and low-rank covariance via iterative application of GLASSO and RCA.
sigma2_n = .01*trace(Cy); % 0.3 * trace(Cy)/d;        % Noise variance. (0.3)
[S D] = eig(Cy);     [D perm] = sort(diag(D),'descend');    % Initialise W with a PCA low-rank estimate.
W_hat_old = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));
WWt_hat_old = W_hat_old * W_hat_old';
Lambda_hat_old = pdinv(Cy);                                 % Initialise Lambda_hat with the empirical inverse-covariance.
%     Lambda_hat_old = eye(d); % *.01*j ;
tic
parfor (i = 1:length(lambda),8)            % Try different magnitudes of lambda.
    [WWt_hat_new, Lambda_hat_new, Lambda_hat_new_inv] = emrca(Y, WWt_hat_old, Lambda_hat_old, sigma2_n, lambda(i), nonZero, limit, emrca_options);
    % Plot results.
    figure(5), clf, colormap('hot')
    subplot(131), imagesc(Lambda_hat_new), colorbar, title([ 'EM/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    subplot(132), imagesc(Lambda_hat_new_inv), colorbar, title('\Sigma_{hat}'), colorbar
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
figure(3), hold on, plot(Recalls, Precisions), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), title('Recall-Precision')
text(Recalls+jit, Precisions+jit, num2cell(lambda))
hold on, plot([1,.87,.27],[.275,.26,.57], 'g-', [1,.67,.2], [.275,.3,.5], 'b-'), legend('EM-RCA','Kronecker-Glasso','Glasso') % literature performance
figure(4), hold on, plot(FPRs, Recalls, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'EM/RCA auc: ' num2str(AUC) ], 4), title('ROC');

% end
