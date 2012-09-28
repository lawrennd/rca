function pstats = emrcaRocStats(Lambda, Lambda_hat, threshold)

% EMRCAROCSTATS Compare a connectivity matrix to the ground truth matrix.
%
% FORMAT
% DESC Takes the estimated and ground truth inverse-covariance and computes
% basic basic statistics required for plotting ROC and precision-recall
% curves.
%
% ARG Lambda : the ground truth connectivity matrix.
%
% ARG Lambda_hat : the estimated connectivity matrix.
%
% ARG threshold : ignores all elements whose absolute value is smaller than
% the given threshold.
%
% RETURN pstats : a row vector with format:
%   [#TruePositives #FalsePositives #FalseNegatives #TrueNegatives]
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

