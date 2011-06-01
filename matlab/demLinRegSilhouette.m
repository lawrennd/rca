% dataSetName = 'silhouette';

% Load data.
% [X Y XTest yTest] = mapLoadData(dataSetName);

% [n sx] = size(X);
% sy = size(Y, 2);

Xo = [X ones(n, 1)];
W = pdinv(Xo'*Xo)*Xo'*Y; % Train W.

XoTest = [XTest ones(size(XTest,1), 1)];
YPred_LR = XoTest*W;

Ysqdiff = (YPred_LR - YTest).^2; % Prediction error.
RMSerror_LR = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff));

% xyzankurAnimCompare(YPred_LR, YTest, 96); % Visualise recovered Y.