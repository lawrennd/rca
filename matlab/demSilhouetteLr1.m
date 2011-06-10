% DEMSILHOUETTELR1 Linear regression demo on Agarwal & Triggs silhouette data.
%
% FORMAT
% DESC Does a least squares regression fit on the weights for Agarwal &
% Triggs silhouette data and animates the result.
%
% SEEALSO : xyzankurAnimCompare
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011
%
% RCA

if ~exist('X','var')
    dataSetName = 'silhouette';
    [X Y XTest yTest] = mapLoadData(dataSetName); % Load data.
    [n sx] = size(X);
    sy = size(Y, 2);
end

Xo = [X ones(n, 1)];
W = pdinv(Xo'*Xo)*Xo'*Y; % Train W.

XoTest = [XTest ones(size(XTest,1), 1)];
YPred_LR = XoTest*W;

Ysqdiff = (YPred_LR - YTest).^2; % Prediction error.
RMSerror_LR = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff));

% xyzankurAnimCompare(YPred_LR, {YTest}, 96, {'YTest'}); % Visualise recovered Y.