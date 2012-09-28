function [W_hat_final, Lambda_hat_final, Lambda_hat_inv_final] = emrca(...
    Y, WWt_hat_old, Lambda_hat_old, sigma2_n, lambda, nonZero, options)

% EMRCA Through a hybrid EM-RCA iterative algorith, learns a Gaussian
% model, whose covariance is composed of a low-rank and a sparse-inverse
% component. The M step is done with the Glasso optimisation algorithm.
% 
% options:
%   verbose: level of verbosity (0: no output, 1: iter (default))
%   showProgress: plot lml function per iteration (0: no plotting, 1: plot lml)
%   errorCheck: check for errors in convergence (0: no error checking, 1:
%   check for errors)
%
% FORMAT
% DESC
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2012
%
% RCA

% subindex = @(x,n,m) x(n,m);

%% Recovery of sparse-inverse and low-rank covariance via iterative
% application of EM (with GLASSO) and RCA.
[n, d] = size(Y);
Cy = Y' * Y /n;
Lambda_hat_new = eye(d);
L1GPoptions.order = -1;         % -1: L-BFGS (limited-memory), 1: BFGS (full-memory), 2: Newton
L1GPoptions.verbose = 0;
    % warmInit = true;
if options.showProgress
    if options.errorCheck
        figure(2), clf
    end
    figure(3), clf
    semilogy(0, options.limit, '.'), xlim([0 options.maxNumIter]), hold on
%     ylim([options.limit subindex(ylim,1,2)])
end
iter = 1;
lml_old = -Inf;
lowerBound_m = -Inf;
lml_new_rca = -Inf;
converged = false;
while ~converged && iter < options.maxNumIter
    if options.verbose
        fprintf('\nEM-RCA iteration: %d\n', iter);
    end
    
    %% E step.
    WWt_plus_noise_inv = pdinv( WWt_hat_old + sigma2_n*eye(d) );
    V_f = pdinv(  WWt_plus_noise_inv  +  Lambda_hat_old ); % Posterior variance of f.
    E_f = (V_f * WWt_plus_noise_inv *  Y')'; % Posterior expectations E[f_n] as rows.
    Avg_E_fft = V_f  +  (E_f'*E_f)./n; % Second moment expectation.
        %             if rank(Avg_E_fft) < d
        %                 warning([ 'rank(Avg_E_fft) = ', num2str(rank(Avg_E_fft)) ]); %#ok<*WNTAG>
        %             end

    if (options.errorCheck)
        %% Variational lower bound after E step. *Should equal the lml*.
            %             [log(det(Lambda_hat_old))*n/2 - sum(sum(Avg_E_fft.*Lambda_hat_old))*n/2 - lambda*sum(abs(Lambda_hat_old(:)))*n/2]
        [lowerBound_e,~,~] = computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat_old, lambda);
        Theta_hat = WWt_hat_old + pdinv(Lambda_hat_old) + sigma2_n*eye(d);
        lml_new_e = computeLogMarginalLikelihood(Cy, n, d, Theta_hat, Lambda_hat_old, lambda);
        
        if (lowerBound_e - lml_new_rca < -1e-6)
            warning([num2str(lowerBound_e - lml_new_rca) ' LowerBound smaller than LML_rca after this E step! Iter:' num2str(iter)]);
            break
        end
        if (abs(lml_new_e - lowerBound_e) > 1e-6) % too high?! was 1e-8
            warning([num2str(lml_new_e - lowerBound_e) ' significant difference between LML and LB_e after this E step! Iter:' num2str(iter)]); %#ok<*WNTAG>
            break
        end
        if (lowerBound_e - lowerBound_m < -1e-6)  % too high?! was 0
            warning([num2str(lowerBound_e - lowerBound_m) ' drop of LowerBound after this E step ! Iter:' num2str(iter)]);
            break
        end
    end
    
    %% M step. Maximise p(f|Lambda) wrt Lambda, via GLASSO.
    warmLambda_hat = Lambda_hat_old;
    funObj = @(x)sparsePrecisionObj(x, d, nonZero, Avg_E_fft);
        %             computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, warmLambda_hat, lambda)
    Lambda_hat_new(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda*ones(d*d,1), L1GPoptions);
        %             computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat_new, lambda(i))
    try
        Lambda_hat_inv_new = pdinv(Lambda_hat_new);
    catch %#ok<CTCH>
        warning(['Matrix is non positive definite tried 10 times adding jitter Iter:' num2str(iter)]);
        break
    end
        %     S = ( Y'*Y - Y'*E_f - E_f'*Y + n*Avg_E_fft ) ./ n;    % Optimise sigma numerically.
        %     funObj_sigma = @(x)Ly_f(x, WWt_hat_old, S, n);
        %     sigma2_n = minConf_TMP(funObj_sigma, sigma2_n, 0, trace(Cy), L1GPoptions);

    %% EM feedback.
    if (options.errorCheck)
        Theta_hat = WWt_hat_old + sigma2_n*eye(d) + Lambda_hat_inv_new;
        lml_new_em = computeLogMarginalLikelihood(Cy, n, d, Theta_hat, Lambda_hat_new, lambda);
        
        % Variational lower bound after M step. ** Should be less than the lml and increased **.
        [lowerBound_m,~,~] = computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat_new, lambda);
        
        if (lml_new_em - lowerBound_m < -1e-6)  % too high?! was -1e-9
            warning([num2str(lml_new_em - lowerBound_e) ' LML smaller than LB_m after this M step! Iter:' num2str(iter)]);
            break
        end
        if (lowerBound_m - lowerBound_e < -1e-6) % !!! too high?
            warning([num2str(lowerBound_m - lowerBound_e) ' LB_m smaller than LB_e after this M step! Iter:' num2str(iter)]);
            break
        end
        
        if options.showProgress
            figure(2),  hold on, plot(  iter, lml_new_em,'.b', ...
                iter-.1, lml_new_e,'.b', ...
                iter-.2, lowerBound_m,'.g', ...
                iter-.3, lowerBound_e,'.g')
        end
    end
    
    %% RCA step: Maximisation of p(y) wrt to W, Cy partially explained by Lambda.
    Theta_explained = Lambda_hat_inv_new + sigma2_n*eye(d);
    [S D] = eig(Cy, Theta_explained);    [D perm] = sort(diag(D),'descend');
    W_hat_new = Theta_explained * S(:,perm(D>1)) * sqrt(diag(D(D>1)-1));
    WWt_hat_new = W_hat_new * W_hat_new';
    
    %% RCA feedback
    Theta_hat = WWt_hat_new + Lambda_hat_inv_new + sigma2_n*eye(d);
    lml_new_rca = computeLogMarginalLikelihood(Cy, n, d, Theta_hat, Lambda_hat_new, lambda);

	% RCA error check.
    if (options.errorCheck)
        if options.showProgress
            figure(2), plot(iter+.5, lml_new_rca,'.r', iter, lml_new_em,'.b'), hold on
        end
        if (lml_new_rca - lml_new_em < -1e-9)
            warning([num2str(lml_new_rca - lml_new_em) ' lml drop observed after RCA iteration! Iter:' num2str(iter)]);
            break;
        end
    end
    
    % Convergence / error check.
    figure(3), semilogy(iter, lml_new_rca - lml_old, 'xr')
    if (lml_new_rca - lml_old) < options.limit
        if (lml_new_rca - lml_old < -1e-6) % too high?!
            warning([num2str(lml_new_rca - lml_old) ' lml drop observed after this iteration! Iter:' num2str(iter)]);
            break
        else
            converged = true;
            fprintf('EM-RCA converged.\n\n')
        end
    end
    
    %% Prepare for new iteration.
    lml_old = lml_new_rca;
        %     warmInit = true;
    Lambda_hat_old = Lambda_hat_new;    WWt_hat_old = WWt_hat_new;
    Lambda_hat_final = Lambda_hat_new;  W_hat_final = W_hat_new;
    Lambda_hat_inv_final = Lambda_hat_inv_new;
    iter = iter + 1;
end

if iter == options.maxNumIter
    fprintf('Maximum number of iterations reached.\n')
	Lambda_hat_final = Lambda_hat_new;
    W_hat_final = W_hat_new;
    Lambda_hat_inv_final = Lambda_hat_inv_new;
end

