function [WWt_hat_new, Lambda_hat_new, Lambda_hat_new_inv] = emrca(Y, WWt_hat_old, Lambda_hat_old, sigma2_n, n, lambda, limit)

% EMRCA Learns an additive Gaussian model, made of a low-rank and a
% sparse component, through a hybrid EM-RCA iterative algorithm. The M step
% is performed via Glasso optimisation.
%
% FORMAT
% DESC
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011, 2012
%
% RCA


%% Recovery of sparse-inverse and low-rank covariance via iterative
% application of EM (with GLASSO) and RCA.
d = size(WWt_hat_old, 1);
Cy = Y' * Y /n;
Lambda_hat_new = eye(d);
nonZero = find(ones(d));    % To induce any prior knowledge of zeros. Typically all ones.
options.order = -1;         % -1: L-BFGS (limited-memory), 1: BFGS (full-memory), 2: Newton
options.verbose = 0;
    % warmInit = true;
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
    [lowerBound_e,~,~] = computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat_old, lambda);
    Theta_hat = WWt_hat_old + sigma2_n*eye(d) + pdinv(Lambda_hat_old);
    lml_new_e = computeLogMarginalLikelihood(Cy, n, d, Theta_hat, Lambda_hat_old, lambda);
    if (abs(lml_new_e - lowerBound_e) > 1e-9)
        warning([num2str(lml_new_e - lowerBound_e) ' significant difference between LML and LB_e after this E step !']); %#ok<*WNTAG>
    end
    if (lowerBound_m > lowerBound_e)
        warning([num2str(lowerBound_e - lowerBound_m) ' drop of the LB after this E step !']);
        break
    end
    
    %% M step. Maximise p(f|Lambda) wrt Lambda, via GLASSO.
    warmLambda_hat = Lambda_hat_old;    % warmSigma_hat = Sigma_hat_old;
        %             [Sigma_hat_new, Lambda_hat_new, iter, avgTol, hasError] = ...
        %                 glasso ( d, Avg_E_fft, 0, lambda(i).*ones(d), ...   % numVars, empirical covariance, computePath, regul.matrix
        %                 0, warmInit, 1, 1, ...  % approximate, warmInit, verbose, penalDiag
        %                 1e-4, 1e2, ...          % tolThreshold (1e-4), maxIter (1e2)
        %                 warmSigma_hat, warmLambda_hat );
    funObj = @(x)sparsePrecisionObj(x, d, nonZero, Avg_E_fft);
        %             computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, warmLambda_hat, lambda)
    Lambda_hat_new(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda*ones(d*d,1), options);
        %             computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat_new, lambda(i))
        %             if any(asym(Lambda_hat_new))
        %             warning([ 'GLasso produced asymmetric Lambda_hat_new by ',...
        %                 num2str(asym(Lambda_hat_new)), '. Lambda_hat_new not symmetrified.' ]);
        %                 Lambda_hat_new = (Lambda_hat_new + Lambda_hat_new') ./ 2; %   Symmetrify.
        %             end
    Lambda_hat_new_inv = pdinv(Lambda_hat_new);
    
    %% Variational lower bound after M step. *Should be less than the lml and increased*.
    [lowerBound_m,~,~] = computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat_new, lambda);
    
    %% EM feedback.
    Theta_hat = WWt_hat_old + sigma2_n*eye(d) + Lambda_hat_new_inv;
    lml_new_em = computeLogMarginalLikelihood(Cy, n, d, Theta_hat, Lambda_hat_new, lambda);
    if (lml_new_em - lowerBound_m < -1e-9)
        warning([num2str(lml_new_em - lowerBound_e) ' LML smaller than LB_m after this M step !']);
        break
    end
    if (lowerBound_e > lowerBound_m)
        warning([num2str(lowerBound_m - lowerBound_e) ' LB_m smaller than LB_e after this M step !']);
        break
    end
        %{
                    fprintf( ['GLasso:\n ' ...
                        ... 'GLasso iterations: %d\n '...
                        ... 'Lambda_hat_new assymetry: %f\n ' ...
                        ... 'avgTol: %e\n hasError: %d\n '...
                        'lambda: %e\n lml_new after EM: %f\n'], ...
                        ... iter, ...
                        ... asym(Lambda_hat_new), ...
                        ... avgTol, hasError, ...
                        lambda(i), lml_new_em );
        %}
    
    %{
    figure(2), plot( k, lml_new_em,'.b', ...
        k-.1, lml_new_e,'.b', ...
        k-.2, lowerBound_m,'.g', ...
        k-.3, lowerBound_e,'.g'), hold on
    %}
    
    %% RCA step: Maximisation of p(y) wrt to W, Cy partially explained by Lambda.
    Theta_explained = Lambda_hat_new_inv + sigma2_n*eye(d);
    [S D] = eig(Cy, Theta_explained);    [D perm] = sort(diag(D),'descend');
    W_hat_new = Theta_explained * S(:,perm(D>1)) * sqrt(diag(D(D>1)-1));
    WWt_hat_new = W_hat_new * W_hat_new';
    
    %% RCA feedback
    Theta_hat = WWt_hat_new + Lambda_hat_new_inv + sigma2_n*eye(d);
    lml_new_rca = computeLogMarginalLikelihood(Cy, n, d, Theta_hat, Lambda_hat_new, lambda);
        %         fprintf('RCA:\n rank(WWt_hat_new): %d\n lml_new after RCA: %f\n\n', ...
        %             rank(WWt_hat_new), lml_new_rca);
        
    %     figure(2), plot(k+.5, lml_new_rca,'.r', k, lml_new_em,'.b'), hold on
    
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
        %     warmInit = true;
    Lambda_hat_old = Lambda_hat_new;    WWt_hat_old = WWt_hat_new;  % Sigma_hat_old = Sigma_hat_new;
    k = k + 1;
    
        %{
        % Plot results of this iteration.
        figure(5), clf, colormap('hot')
        subplot(131), imagesc(Lambda_hat_new), colorbar
        title([ 'GLasso/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
        subplot(132), imagesc(Lambda_hat_new_inv), colorbar, title('\Sigma_{hat}'), colorbar
        subplot(133), imagesc(WWt_hat_new), colorbar, title('RCA-recovered WW'''), colorbar
        %}
end

