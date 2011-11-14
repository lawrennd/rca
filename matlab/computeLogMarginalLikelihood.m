function lml = computeLogMarginalLikelihood(C, n, d, Sigma)

try
    [Sigma_inv, R] = pdinv(Sigma);
    % Fast Way to compute -logdet(C) - tr(C*Sigma)
    lml = -log(2*pi)*d*n/2 - 2*sum(log(diag(R)))*n/2 - sum(sum(n*C .* Sigma_inv))/2;
    %     -log(2*pi)*d*n/2 - log(det(Sigma))*n/2 - sum(sum(n*C.*pdinv(Sigma)))/2
catch %#ok<*CTCH>
    error('Non positive definite covariance.')
end