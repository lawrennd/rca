% function [lb, Q, H] = computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat)
function [lb, Q, H] = computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat)

[n d] = size(Y);

Cyy = Y' * Y / n;
Cfy = E_f' * Y / n;
Cff = E_f' * E_f / n;

WWt_s2 = WWt_hat_old + sigma2_n*eye(d);
[WWt_s2_inv, WWt_s2_CU] = pdinv( WWt_s2 );
[~, Lambda_hat_CU] = pdinv( Lambda_hat );
[V_f, V_f_inv_CU] = pdinv( WWt_s2_inv + Lambda_hat_old );


% Expectation of complete-data log-likelihood, wrt posterior f|y.
Q = - log(2 * pi) * d * n ...
    - log(det(WWt_s2)) * n/2 ...
    + log(det(Lambda_hat)) * n/2 ...
    - sum(sum( (Cyy - 2*Cfy)' .* WWt_s2_inv )) * n/2 ...
    - sum(sum( (pdinv(WWt_s2_inv + Lambda_hat_old) + Cff)' .* (WWt_s2_inv + Lambda_hat) ))  * n/2;

    - log(2 * pi) * d * n ...
    - 2*sum(log(diag(WWt_s2_CU))) * n/2 ...
    + 2*sum(log(diag(Lambda_hat_CU))) * n/2 ...
    - sum(sum( (Cyy - 2*Cfy)' .* WWt_s2_inv )) * n/2 ...
    - sum(sum( (V_f + Cff) .* (WWt_s2_inv + Lambda_hat) ))  * n/2

% Entropy of posterior f|y.
H = ( log(2 * pi) * d - log(det(Lambda_hat_old + WWt_s2_inv)) + d ) * n/2;
    ( log(2 * pi) * d - 2*sum(log(diag(V_f_inv_CU))) + d ) * n/2

lb = Q + H;