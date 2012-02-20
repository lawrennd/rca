function pstats = emrcaRocStats(Lambda, Lambda_hat, threshold)

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

if nargin < 3 || isempty(threshold)
    threshold = 0;
end

%% Performance stats.
A = boolean(triu(Lambda,1) ~= 0);               % Binary matrix indicating edges in the ground truth Lambda.
triuLambda_hat = triu(Lambda_hat, 1);
B = boolean( abs(triuLambda_hat) > threshold);  % Edges in the estimated Lambda.

%% Row format: [ TP FP FN TN ]
TP = sum( A(:) & B(:) );
FP = sum( ~A(:) & B(:) );
FN = sum( A(:) & ~B(:) );
TN = sum( ~A(:) & ~B(:) );
pstats = [ TP FP FN TN ];

