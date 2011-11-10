function [lb, Q, H] = computeLowerBound(Y, E_f, WWt_hat_old, Lambda_hat_old, sigma2_n, Lambda_hat)

Cy = Y'*Y / n;
Cfy = E_f'*Y / n;
Cf = E_f'*E_f / n;

WWt_s2 = WWt_hat_old + sigma2_n*eye(d);
WWt_s2_inv = pdinv(WWt_s2);

% Expectation of complete-data log-likelihood, wrt posterior f|y.
Q = - log(2 * pi) * d * n ...
    - log(det(WWt_s2_inv * Lambda_hat) ) * n/2 ...
    - ( trace((Cy - 2*Cfy + Cf) * WWt_s2_inv) ...
        + trace(pdinv(WWt_s2_inv + Lambda_hat_old) * (WWt_s2_inv + Lambda_hat)) ...
        + trace(Cf * Lambda_hat) ) * n/2;

% Entropy of posterior f|y.
H = ( log(2 * pi) * d - log(det(Lambda_hat_old + WWt_s2_inv)) + n ) * n/2;

lb = Q + H;