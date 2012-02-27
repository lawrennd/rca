% GRAD_LY_SIGMA2_N Gradient of expected complete data log-likelihood wrt
% sigma2_n.
%
% FORMAT
% DESC Gradient of expected log N(Y|F,WWt + sigma2_n*I) wrt
% sigma2_n. Expectation is wrt the posterior p(F|Y).
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011
%
% RCA

function gradient = Grad_Ly_f_sigma2_n(sigma, WWt, S, n)

d = size(WWt,1);

% gradient = zeros(size(sigma));
% if numel(sigma) > 1
%     for i = 1:size(sigma,1)
%         for j = 1:size(sigma,2)
%             Sigma_inv = pdinv( WWt + sigma(i,j)*eye(d) );
%             gradient(i,j) = -trace(Sigma_inv  -  Sigma_inv * S * Sigma_inv) * n/2; % Gradient of Ly_f wrt to sigma2_n.
%         end
%     end
% else
    Sigma_inv = pdinv( WWt + sigma*eye(d) );
    gradient = trace(Sigma_inv  -  Sigma_inv * S * Sigma_inv) * n/2; % Gradient of Ly_f wrt to sigma2_n.
% end

