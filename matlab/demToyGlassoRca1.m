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

clear, clc
addpath(genpath('~/mlprojects/matlab/general/'))
addpath(genpath('~/mlprojects/rca/matlab/glasso/'))
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
importTool({'rca','ndlutil','gprege'})
asym = @(x) sum(sum(abs(x - x.'))); % Asymmetry test.

figure(1), clf, colormap('hot')
figure(2), clf, colormap('hot')
figure(3), clf, colormap('hot')
figure(4), clf, colormap('hot')
figure(5), clf, colormap('hot')

limit = 1e-4;
lambda = 5.^linspace(-8,3,30);
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
%     testLambda = Lambda + diag(abs(randn(d,1))).*4; % 4
    [~, pd] = chol(testLambda);
    diagonal = diagonal*2;
end
Lambda = testLambda;    Sigma = pdinv(Lambda);
% Sample from p(y).
W = randn(d,p);
WWt = W*W';
Theta = WWt + Sigma + sigma2_n*eye(d);
Y = gaussSamp(Theta, n);
figure(1), clf, colormap('hot')
subplot(131), imagesc(Lambda), title('sparse \Lambda'), colorbar
subplot(132), imagesc(Sigma), title('sparse-inverse \Sigma'), colorbar
subplot(133), imagesc(WWt), title('Low-rank WW'''), colorbar

Y = Y - repmat(mean(Y),n,1);
Cy = Y' * Y /n;
%}


%% Standard Glasso on (un)confounded simulated data, with varying lambda.
%{
confounders{1} = WWt;   confounders{2} = zeros(size(WWt));
linestyle = {'-xb','--om'}; legends = {'GL','"Ideal" GL'};
AUCs = zeros(2,1); figure(2), clf, figure(4), clf
for c = 1:2
    s = RandStream('mcg16807','Seed', 666); RandStream.setDefaultStream(s) % 23, 1e5
    Y_ = gaussSamp(confounders{c} + Sigma + sigma2_n*eye(d), n);   Y_ = Y_ - repmat(mean(Y_),n,1);
    Cy_ = Y_' * Y_/n;
	nonZero = find(ones(d));    % To induce any prior knowledge of non-zeros. Typically all ones.
    warmLambda_hat = Lambda; % eye(d);
    funObj = @(x)sparsePrecisionObj(x, d, nonZero, Cy_);
    options.verbose = 1;    options.order = -1;
    A = boolean( triu(Lambda,1) ~= 0 );
    for i = 1:length(lambda)
        Lambda_hat{i} = eye(d);
        Lambda_hat{i}(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda(i)*ones(d*d,1), options);
        %         [Sigma_hat{i}, Lambda_hat{i}] = ...
        %                 glasso( d, Cy_, 0, lambda(i)*ones(d), ...   % numVars, empirical covariance, computePath, regul.matrix
        %                 0, 0, 0, 1, ...                             % approximate, warmInit, verbose, penalDiag
        %                 1e-4, 1e4, ...                              % tolThreshold (1e-4), maxIter (1e4)
        %                 zeros(d), zeros(d));                        % warmLambda, warmSigma
        triuLambda_hat{i} = triu(Lambda_hat{i}, 1);
        figure(3), imagesc(Lambda_hat{i}), colormap(hot), colorbar,...
            title([ 'GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
        % Evaluation
        B{i} = boolean( triuLambda_hat{i} ~= 0 );
        TPs(i) = sum( A(:) & B{i}(:) );
        FPs(i) = sum( ~A(:) & B{i}(:) );
        FNs(i) = sum( A(:) & ~B{i}(:) );
        TNs(i) = sum( ~A(:) & ~B{i}(:) );
    end
    Recalls = TPs ./ (TPs + FNs);
    Precisions = TPs ./ (TPs + FPs);
    FPRs = FPs ./ (FPs + TNs);
    figure(2), hold on, plot(Recalls, Precisions, linestyle{c}), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
    figure(4), hold on, plot(FPRs, Recalls, linestyle{c}), xlim([0 1]), xlabel('FPR'), ylabel('TPR')
    AUCs(c) = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);
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

% lambda = 10^-2;           % 10^-2 too strong
for i = 1:length(lambda)    % Try different magnitudes of lambda.
    
    % Initialise W with a PPCA low-rank estimate.
    [S D] = eig(Cy);     [D perm] = sort(diag(D),'descend');
    W_hat_old = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));

%     W_hat_old = zeros(d,p);
%     W_hat_old = W;  % Use the true low-rank.
    
    WWt_hat_old = W_hat_old * W_hat_old';
    
%     Lambda_hat_old = Lambda;    Sigma_hat_old = pdinv(Lambda);
    Lambda_hat_old = eye(d);    Sigma_hat_old = pdinv(Lambda);
    Lambda_hat_new = eye(d);    Sigma_hat_new = eye(d);
    
    nonZero = find(ones(d));    % To induce any prior knowledge of zeros. Typically all ones.
    options.order = -1;         % -1: L-BFGS (limited-memory), 1: BFGS (full-memory), 2: Newton
    options.verbose = 0;
    warmInit = true;
    figure(2), clf
    k = 1;
    lml_old = -Inf;
    lowerBound_m = -Inf;
    converged = false;
    while ~converged
        fprintf('\nEM-RCA iteration: %d\n', k);

        %% E step.
        WWt_plus_noise_inv = pdinv( WWt_hat_old + sigma2_n*eye(d) );
        V_f = pdinv(  WWt_plus_noise_inv  +  Lambda_hat_old ); % Posterior variance of f.
        E_f = (V_f * WWt_plus_noise_inv *  Y')'; % Posterior expectations E[f_n] as rows.
        Avg_E_fft = V_f  +  (E_f'*E_f)./n; % Second moment expectation.
            %             if rank(Avg_E_fft) < d
            %                 warning([ 'rank(Avg_E_fft) = ', num2str(rank(Avg_E_fft)) ]); %#ok<*WNTAG>
            %             end

        %% Variational lower bound after E step. *Should equal the lml*.
            %             [log(det(Lambda_hat_old))*n/2 - sum(sum(Avg_E_fft.*Lambda_hat_old))*n/2 - lambda*sum(abs(Lambda_hat_old(:)))*n/2]
        [lowerBound_e, Q_e, H_e] = computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat_old, lambda(i));
        Theta_hat = WWt_hat_old + sigma2_n*eye(d) + pdinv(Lambda_hat_old);
        lml_new_e = computeLogMarginalLikelihood(Cy, n, d, Theta_hat, Lambda_hat_old, lambda(i));
        %             lml_new_e = -log(2*pi)*d*n/2 - log(det(Theta_hat))*n/2 - sum(sum((Y'*Y)'.*pdinv(Theta_hat)))/2;
        if (abs(lml_new_e - lowerBound_e) > 1e-9)
            warning([num2str(lml_new_e - lowerBound_e) ' significant difference between LML and LB_e after this E step !']); %#ok<*WNTAG>
        end
        if (lowerBound_m > lowerBound_e)
            warning([num2str(lowerBound_e - lowerBound_m) ' drop of the LB after this E step !']);
            break
        end

        %% M step. Maximise p(f|Lambda) wrt Lambda, via GLASSO.
        warmLambda_hat = Lambda_hat_old;    warmSigma_hat = Sigma_hat_old;
            %             [Sigma_hat_new, Lambda_hat_new, iter, avgTol, hasError] = ...
            %                 glasso ( d, Avg_E_fft, 0, lambda(i).*ones(d), ...   % numVars, empirical covariance, computePath, regul.matrix
            %                 0, warmInit, 1, 1, ...  % approximate, warmInit, verbose, penalDiag
            %                 1e-4, 1e2, ...          % tolThreshold (1e-4), maxIter (1e2)
            %                 warmSigma_hat, warmLambda_hat );
        funObj = @(x)sparsePrecisionObj(x, d, nonZero, Avg_E_fft);
            %             [log(det(warmLambda_hat))*n/2 - sum(sum(Avg_E_fft.*warmLambda_hat))*n/2 - lambda*sum(abs(warmLambda_hat(:)))*n/2]
            %             computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, warmLambda_hat, lambda)
        Lambda_hat_new(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda(i)*ones(d*d,1), options);
            %             [log(det(Lambda_hat_new))*n/2 - sum(sum(Avg_E_fft.*Lambda_hat_new))*n/2 - lambda(i)*sum(abs(Lambda_hat_new(:)))*n/2]
            %             computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat_new, lambda(i))
            %             if any(asym(Lambda_hat_new))
            %             warning([ 'GLasso produced asymmetric Lambda_hat_new by ',...
            %                 num2str(asym(Lambda_hat_new)), '. Lambda_hat_new not symmetrified.' ]);
            %                 Lambda_hat_new = (Lambda_hat_new + Lambda_hat_new') ./ 2; %   Symmetrify.
            %             end
        Lambda_hat_new_inv = pdinv(Lambda_hat_new);

        %% Variational lower bound after M step. *Should be less than the lml and increased*.
        [lowerBound_m, Q_m, H_m] = computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat_new, lambda(i));

        %% EM feedback.
        Theta_hat = WWt_hat_old + sigma2_n*eye(d) + Lambda_hat_new_inv;
        lml_new_em = computeLogMarginalLikelihood(Cy, n, d, Theta_hat, Lambda_hat_new, lambda(i));
        %             lml_new_em = -log(2*pi)*d*n/2 - log(det(Theta_hat))*n/2 - sum(sum((Y'*Y)'.*pdinv(Theta_hat)))/2;
        if (lml_new_em - lowerBound_m < -1e-9)
            warning([num2str(lml_new_em - lowerBound_e) ' LML smaller than LB_m after this M step !']);
            break
        end
        if (lowerBound_e > lowerBound_m)
            warning([num2str(lowerBound_m - lowerBound_e) ' LB_m smaller than LB_e after this M step !']);
            break
        end

            %             fprintf( ['GLasso:\n ' ...
            %                 ... 'GLasso iterations: %d\n '...
            %                 ... 'Lambda_hat_new assymetry: %f\n ' ...
            %                 ... 'avgTol: %e\n hasError: %d\n '...
            %                 'lambda: %e\n lml_new after EM: %f\n'], ...
            %                 ... iter, ...
            %                 ... asym(Lambda_hat_new), ...
            %                 ... avgTol, hasError, ...
            %                 lambda(i), lml_new_em );
        figure(2), plot( k, lml_new_em,'.b', ...
            k-.1, lml_new_e,'.b', ...
            k-.2, lowerBound_m,'.g', ...
            k-.3, lowerBound_e,'.g'), hold on

        %% RCA step: Maximisation of p(y) wrt to W, Cy partially explained by Lambda.
        Theta_explained = Lambda_hat_new_inv + sigma2_n*eye(d);
        [S D] = eig(Cy, Theta_explained);    [D perm] = sort(diag(D),'descend');
        W_hat_new = Theta_explained * S(:,perm(D>1)) * sqrt(diag(D(D>1)-1));
        WWt_hat_new = W_hat_new * W_hat_new';
        
        %% RCA feedback
        Theta_hat = WWt_hat_new + Lambda_hat_new_inv + sigma2_n*eye(d);
        lml_new_rca = computeLogMarginalLikelihood(Cy, n, d, Theta_hat, Lambda_hat_new, lambda(i));
%         -log(2*pi)*d*n/2 - log(det(Theta_hat))*n/2 - sum(sum((Y'*Y)'.*pdinv(Theta_hat)))/2 - (lambda(i) * sum(abs(Lambda_hat_new(:))))*n/2
%         fprintf('RCA:\n rank(WWt_hat_new): %d\n lml_new after RCA: %f\n\n', ...
%             rank(WWt_hat_new), lml_new_rca);
        figure(2), plot(k+.5, lml_new_rca,'.r', k, lml_new_em,'.b'), hold on
        
        % RCA error check.
        if lml_new_rca < lml_new_em
            warning([num2str(lml_new_rca - lml_new_em) ' lml drop observed after RCA iteration!']);
            break;
        end
        
        % Convergence / error check.
        if (lml_new_rca - lml_old) < limit
            if lml_old > lml_new_rca
                warning([num2str(lml_new_rca - lml_old) ' lml drop observed after this iteration!']);
                break
            else
                converged = true;
                fprintf('EM-RCA algorithm converged.\n\n')
            end
        end
        
        %% Prepare for new iteration.
        lml_old = lml_new_rca;
        warmInit = true;
        Lambda_hat_old = Lambda_hat_new;    WWt_hat_old = WWt_hat_new;  Sigma_hat_old = Sigma_hat_new;
        k = k + 1;
        
            % Plot results of this iteration.
            %         figure(5), clf, colormap('hot')
            %         subplot(131), imagesc(Lambda_hat_new), colorbar
            %             title([ 'GLasso/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
            %         subplot(132), imagesc(Lambda_hat_new_inv), colorbar, title('\Sigma_{hat}'), colorbar
            %         subplot(133), imagesc(WWt_hat_new), colorbar, title('RCA-recovered WW'''), colorbar
    end
    
    % Plot results.
    figure(5), clf, colormap('hot')
    subplot(131), imagesc(Lambda_hat_new), colorbar
        title([ 'GLasso/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    subplot(132), imagesc(Lambda_hat_new_inv), colorbar, title('\Sigma_{hat}'), colorbar
    WWt_hat_new(WWt_hat_new > max(WWt(:))) = max(WWt(:));   WWt_hat_new(WWt_hat_new < min(WWt(:))) = min(WWt(:));
    subplot(133), imagesc(WWt_hat_new), colorbar, title('RCA-recovered WW'''), colorbar
    
    % Performance stats.
    A = boolean(triu(Lambda,1) ~= 0);
    triuLambda_hat{i} = triu(Lambda_hat_new, 1);
    B{i} = boolean( triuLambda_hat{i} ~= 0 );
    TPs(i) = sum( A(:) & B{i}(:) );
    FPs(i) = sum( ~A(:) & B{i}(:) );
    FNs(i) = sum( A(:) & ~B{i}(:) );
    TNs(i) = sum( ~A(:) & ~B{i}(:) );
end
Recalls = TPs ./ (TPs + FNs);
Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);
AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);
figure(3), hold on, plot(Recalls, Precisions, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), title('Recall-Precision')
figure(4), hold on, plot(FPRs, Recalls, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('FPR'), ylabel('TPR')
legend([ 'RCA-GLasso auc: ' num2str(AUC) ], 4), title('ROC');
%}

