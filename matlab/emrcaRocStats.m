function pstats = emrcaRocStats(Lambda, Lambda_hat_new)

% PERFORMANCESTATS Takes the estimated and ground truth inverse-covariance
% and computes basic statistics required for plotting ROC and
% precision-recall curves.
%
% FORMAT
% DESC
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011, 2012
%
% RCA

%% Performance stats.
A = boolean(triu(Lambda,1) ~= 0);   % Binary matrix indicating dges in the ground truth Lambda.
triuLambda_hat = triu(Lambda_hat_new, 1);
B = boolean( triuLambda_hat ~= 0 );   % Edges in the estimated Lambda.

%% Row format: [ TP FP FN TN ]
TP = sum( A(:) & B(:) );
FP = sum( ~A(:) & B(:) );
FN = sum( A(:) & ~B(:) );
TN = sum( ~A(:) & ~B(:) );
pstats = [ TP FP FN TN ];
