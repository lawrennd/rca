% DEMTOYGLASSORCA1 RCA-Glasso demo on simulated data.
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
addpath(genpath('~/mlprojects/rca/matlab/glasso/'))
importTool({'rca','ndlutil','gprege'})
asym = @(x) sum(sum(abs(x - x.'))); % Asymmetry test.

figure(1), clf, colormap('hot')
figure(2), clf, colormap('hot')
figure(3), clf, colormap('hot')
figure(4), clf, colormap('hot')
figure(5), clf, colormap('hot')

limit = 1e-5;
lambda = 5.^linspace(-8,3,20);
Sigma_hat = cell(length(lambda),1);
Lambda_hat = cell(length(lambda),1);
triuLambda_hat = cell(length(lambda),1);
B = cell(length(lambda),1);
TPs = zeros(length(lambda), 1);
FPs = zeros(length(lambda), 1);
FNs = zeros(length(lambda), 1);
TNs = zeros(length(lambda), 1);


%% Data generation.
d = 50; % Observed dimensions.
p = 3; % Low-rank
n = 100;
sigma2_n = 1e-2; % Noise variance.
s = RandStream('mcg16807','Seed', 66); RandStream.setDefaultStream(s)
validLambda = false;
while ~validLambda
    Lambda = triu(rand(d)<.01, 1) .* (randn(d)*sqrt(2) + 1);
    Lambda = Lambda + Lambda' + diag(abs(randn(d,1))).*4;
    validLambda = all(eig(Lambda)>0);
end
triuLambda = triu(Lambda,1);
Sigma = pdinv(Lambda);
W = randn(d,p);
WWt = W*W';
Theta = WWt + Sigma + sigma2_n*eye(d);
Y = gaussSamp(Theta, n); % Sample from p(y).
figure(1), clf, colormap('hot')
subplot(131), imagesc(Lambda), title('Sparse \Lambda'), colorbar
subplot(132), imagesc(Sigma), title('Sparse-inverse \Sigma'), colorbar
subplot(133), imagesc(WWt), title('Low-rank WW'''), colorbar

Y = Y - repmat(mean(Y),n,1);
Cy = Y'*Y/n; % Sample covariance Y.
totalvar = trace(Cy)/d;


%% Standard Glasso on (un)confounded simulated data, with varying lambda.
%{
confounders{1} = WWt;   confounders{2} = zeros(size(WWt));
linestyle = {'-xb','--om'}; legends = {'GLasso','Ideal GLasso'};
AUCs = zeros(2,1); figure(2), clf, figure(4), clf
for c = 1:2
    s = RandStream('mcg16807','Seed', 666); RandStream.setDefaultStream(s) % 23, 1e5
    Y_ = gaussSamp(confounders{c} + Sigma + sigma2_n*eye(d), n);   Y_ = Y_ - repmat(mean(Y_),n,1);
    Cy_ = Y_'*Y_/n;
    A = boolean(triu(full(Lambda),1)~=0);
    for i = 1:length(lambda)
        [Sigma_hat{i}, Lambda_hat{i}] = glasso(d, Cy_, 0, lambda(i)*ones(d), 0,0,0,0, 1e-4, 1e4, zeros(d), zeros(d));
        triuLambda_hat{i} = triu(Lambda_hat{i}, 1);
        figure(3), imagesc(Lambda_hat{i}), colormap(hot), colorbar,...
            title([ '(RCA)GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
        % Evaluation
        B{i} = boolean( triuLambda_hat{i} ~= 0 );
        TPs(i) = sum( A(:) & B{i}(:) );
        FPs(i) = sum( ~A(:) & B{i}(:) );
        FNs(i) = sum( A(:) & ~B{i}(:) );
        TNs(i) = sum( ~A(:) & ~B{i}(:) );
    end
    TPRs = TPs ./ (TPs + FNs);
    Precisions = TPs ./ (TPs + FPs);
    FPRs = FPs ./ (FPs + TNs);
    figure(2), hold on, plot(TPRs, Precisions, linestyle{c}), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
    figure(4), hold on, plot(FPRs, TPRs, linestyle{c}), xlim([0 1]), xlabel('FPR'), ylabel('TPR')
    AUCs(c) = trapz(flipud(FPRs), flipud(TPRs)) / max(FPRs);
end
figure(2), legend(legends), title('Recall-Precision');
figure(4), legend([ legends{1} ' auc: ' num2str(AUCs(1)) ], [ legends{2} ' auc: ' num2str(AUCs(2)) ], 4), title('ROC');
%}


%% **USELESS** Standard GLasso on confounded simulated data, with varying lambda,
% *after* having explained away the true low-rank structure via RCA.
%{
s = RandStream('mcg16807','Seed', 666); RandStream.setDefaultStream(s) % 23, 1e5
Y_ = gaussSamp(Theta, n);   Y_ = Y_ - repmat(mean(Y_),n,1);
Cy_ = Y_'*Y_/n;
% Retrieve residual variance basis via RCA.
Theta_exp = WWt + sigma2_n*eye(d)*.95;
[S D] = eig(Cy_, Theta_exp);    [D perm] = sort(diag(D),'descend');
V_hat = Theta_exp * S(:,perm(D>1)) * sqrt(diag(D(D>1)-1));
VVt_hat = V_hat * V_hat' + sigma2_n*eye(d)*.05;
figure(3), imagesc(VVt_hat), colorbar
A = boolean(triu(full(Lambda),1) ~= 0);
for i = 1:length(lambda)
    [Sigma_hat{i}, Lambda_hat{i}] = glasso(d, VVt_hat, 0, lambda(i)*ones(d),0,0,0,0,1e-4,1e+4, zeros(d), zeros(d));
    triuLambda_hat{i} = triu(Lambda_hat{i}, 1);
    figure(3), imagesc(Lambda_hat{i}), colormap(hot), colorbar, ...
        title([ '(RCA)GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    % Evaluation
    B{i} = boolean( triuLambda_hat{i} ~= 0 );
    TPs(i) = sum( A(:) & B{i}(:) );
    FPs(i) = sum( ~A(:) & B{i}(:) );
    FNs(i) = sum( A(:) & ~B{i}(:) );
    TNs(i) = sum( ~A(:) & ~B{i}(:) );
end
TPRs = TPs ./ (TPs + FNs);
Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);
figure(2), hold on, plot(TPRs, Precisions, ':rs'), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
title('Recall-Precision');
figure(4), hold on, plot(FPRs, TPRs, ':rs'), xlim([0 1]), xlabel('FPR'), ylabel('TPR')
AUC = trapz(flipud(FPRs), flipud(TPRs)) / max(FPRs);
%}


%% Recovery of low-rank component WWt by explaining away the true
% sparse-inverse covariance Sigma.
%{
Theta_exp = Sigma + sigma2_n*eye(d);
[S D] = eig(Cy, Theta_exp);    [D perm] = sort(diag(D),'descend');
W_hat = Theta_exp * S(:,perm(D>1)) * sqrt(diag(D(D>1)-1));
WWt_hat = W_hat * W_hat';
figure(3), imagesc(WWt_hat), colorbar;
%}


%% Recovery of sparse-inverse and low-rank covariance via iterative
% application of GLASSO and RCA.

lambda = 10^-2;
for i = 1:length(lambda) % Try different magnitudes of lambda.
    % Initialise W with a PPCA low-rank estimate.
    [S D] = eig(Cy);     [D perm] = sort(diag(D),'descend');
    W_hat = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));
    % No prior knowledge on the low-rank component.
    %         W_hat = zeros(d,p);
    %True low-rank.
    W_hat = W;
    WWt_hat = W_hat*W_hat';
    Lambda_hat = zeros(d);
    converged = false;
    lml_old = -Inf; % Initial log-marginal likelihood.
    warmInit = false; warmSigma_hat = zeros(d); warmLambda_hat = zeros(d);
    k = 1;
    figure(2), clf
    while ~converged
        % E step.
        WWt_plus_noise_inv = pdinv( WWt_hat + sigma2_n*eye(d) );
        V_f = pdinv(  WWt_plus_noise_inv  +  Lambda_hat ); % Posterior variance of f.
        E_f = (V_f * WWt_plus_noise_inv *  Y')'; % Posterior expectations E[f_n] as rows.
        Avg_E_fft = V_f  +  (E_f'*E_f)./n; % Second moment expectation.
        if rank(Avg_E_fft) < size(Avg_E_fft,1)
            warning([ 'rank(Avg_E_fft) = ', num2str(rank(Avg_E_fft)) ]);
        end
%{        
        % M step. Maximise E_f|y[ N(Y|F, WW'+sigma2_n*I) ] wrt sigma2_n via grid-search.
        S_f = Cy - (Y'*E_f + E_f'*Y)/n + Avg_E_fft;
        sigmas = logspace(log10(1e-4), log10(totalvar), 100);
        history_L = Ly_f(sigmas, WWt_hat, S_f, n);
        history_GradL = Grad_Ly_f_sigma2_n(sigmas, WWt_hat, S_f, n);
        [~, imax] = max(history_L);
        sigma2_n = sigmas(imax); % !!!
        fprintf('sigma max: %f\ngradient at max: %f\n\n', sigmas(imax), history_GradL(imax) )
%}        
        % M step. Maximise p(f|Lambda) wrt Lambda, via GLASSO.
        [Sigma_hat, Lambda_hat, iter, avgTol, hasError] = glasso(d, Avg_E_fft, 0, lambda(i).*ones(d),...
            0,warmInit,0,1, 1e-4, 1e2, warmSigma_hat, warmLambda_hat);
        Lambda_hat = (Lambda_hat + Lambda_hat') ./ 2; % Symmetrify.
        fprintf('GLasso:\n EM iteration: %d\n GLasso iterations: %d\n Lambda assymetry: %f\n hasError: %d\n\n',...
            k, iter, asym(Lambda_hat), hasError);
        
        Theta_hat = WWt_hat + pdinv(Lambda_hat) + sigma2_n*eye(d);
        lml = -log(2*pi)*d*n/2 - log(det(Theta_hat))*n/2 - sum(sum((Y'*Y)'.*pdinv(Theta_hat)))/2;
        figure(2), plot(k, lml,'.b'), hold on
        
        % RCA step: Maximisation of p(y) wrt to W, Cy partially explained by Lambda.
        Theta_explained = pdinv(Lambda_hat) + sigma2_n*eye(d);                     % Is sigma necessary?
        [S D] = eig(Cy, Theta_explained);    [D perm] = sort(diag(D),'descend');
        W_hat = Theta_explained * S(:,perm(D>1)) * sqrt(diag(D(D>1)-1));
        WWt_hat = W_hat*W_hat';
        fprintf('RCA:\n rank(WWt_hat): %d\n\n', rank(WWt_hat));
        Theta_hat = WWt_hat + pdinv(Lambda_hat) + sigma2_n*eye(d);                 % Is sigma necessary?
        
        lml = -log(2*pi)*d*n/2 - log(det(Theta_hat))*n/2 - sum(sum((Y'*Y)'.*pdinv(Theta_hat)))/2;
        converged = (lml - lml_old) < limit;
        if lml_old > lml
            warning([num2str(lml - lml_old) ' lml decrease']);
        end
        lml_old = lml;
        warmInit = true;    warmSigma_hat = Sigma_hat;  warmLambda_hat = Lambda_hat;
        
        figure(2), plot(k+.01, lml,'.r'), hold on
        k = k + 1;
    end
    % Plot results.
    figure(5), clf, colormap('hot')
    subplot(131), imagesc(Lambda_hat), colorbar, title([ 'GLasso/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    subplot(132), imagesc(Sigma_hat), colorbar, title('\Sigma_{hat}'), colorbar
    WWt_hat(WWt_hat > max(WWt(:))) = max(WWt(:));    WWt_hat(WWt_hat < min(WWt(:))) = min(WWt(:));
    subplot(133), imagesc(WWt_hat), colorbar, title('RCA-recovered WW'''), ylabel('+','fontsize',15), colorbar
    
    % Performance stats.
    triuLambda_hat{i} = triu(Lambda_hat, 1);
    B{i} = boolean( triuLambda_hat{i} ~= 0 );
    TPs(i) = sum( A(:) & B{i}(:) );
    FPs(i) = sum( ~A(:) & B{i}(:) );
    FNs(i) = sum( A(:) & ~B{i}(:) );
    TNs(i) = sum( ~A(:) & ~B{i}(:) );
end
TPRs = TPs ./ (TPs + FNs);
Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);
AUC = trapz(flipud(FPRs), flipud(TPRs)) / max(FPRs);
figure(3), hold on, plot(TPRs, Precisions, '-rs'), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), title('Recall-Precision')
figure(4), hold on, plot(FPRs, TPRs, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('FPR'), ylabel('TPR')
legend([ 'RCA-GLasso auc: ' num2str(AUCs(1)) ], 4), title('ROC');
%}

