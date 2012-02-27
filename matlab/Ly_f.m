% LY_F Negative of the expected complete data log-likelihood.
%
% FORMAT
% DESC The negative expected log N(Y|F,WWt + sigma2_n*I). Expectation is wrt to the
% posterior p(F|Y).
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011, 2012
%
% RCA

function [L, g] = Ly_f(sigma, WWt, S, n)
%% E_f|y[ N(Y|F,WW'+x*I) ]
d = size(WWt,1);

% L = zeros(size(sigma));
% if numel(sigma) > 1
%     for i = 1:size(sigma,1)
%         for j = 1:size(sigma,2)
%             L(i,j) = - (  logdet(WWt + sigma(i,j)*eye(d))  +  trace( S * pdinv(WWt + sigma(i,j)*eye(d)) )  ) * n/2;
%         end
%     end
% else
	L =  (  logdet(WWt + sigma*eye(d))  +  trace( S * pdinv(WWt + sigma*eye(d)) )  ) * n/2;
    g = Grad_Ly_f_sigma2_n(sigma, WWt, S, n);
% end
