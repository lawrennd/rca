% DEMSILHOUETTERCA1 RCA demo on Agarwal & Triggs silhouette data.
%
% FORMAT
% DESC Fits private and shared latent spaces with iterative RCA. Compares
% against linear regression and PCCA.
%
% SEEALSO : pcca, cca, demSilhouetteLr1, xyzankurAnimCompareMultiple
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011
%
% RCA

addpath(genpath('~/mlprojects/matlab/general/'))
importTool({'rca','cca','datasets','ndlutil','mocap'})
fontName = 'times';
fontSize = 26;

% savedState = struct(RandStream.getDefaultStream).State;
% defaultStream.State = savedState;

clf, clc

dataSetName = 'silhouette';
converged = false;
tracesC = [];
k = 1;  t = 1; i = 1;
limit = 1e-3;

% Load data.
[X, Y, XTest, YTest] = mapLoadData(dataSetName);
[n, sx] = size(X);
sy = size(Y, 2);

% Induce isotropic noise.
noise = .03;
X = X + randn(size(X))*noise;
XTest = XTest + randn(size(XTest))*noise;
Y = Y + randn(size(Y))*noise;
YTest = YTest + randn(size(YTest))*noise;

Sy = cov(Y, 1); Sx = cov(X, 1); S = cov([X Y], 1);  % Data covariances.
Sy = Sy + abs(mean(diag(Sy)))*1e-6*eye(sy); S(sx+1:end,sx+1:end) = Sy;  % Add jitter.
Sx = Sx + abs(mean(diag(Sx)))*1e-6*eye(sx); S(1:sx,1:sx) = Sx;

for alpha = .05:.05:1;  % Fraction of mean trace as noise level.
%     alpha = .3; beta = .3;
    noisey = alpha*trace(Sy)/sy;    noisex = alpha*trace(Sx)/sx;
    Ny = noisey*eye(sy); Nx = noisex*eye(sx); Noise = blkdiag(Nx,Ny);  % Fix noise with alpha.
    W1W1 = zeros(sy);   W4W4 = zeros(sx);   % No independent components.
    W2W2 = zeros(sy);   W3W3 = zeros(sx);   % No shared components.
    W2 = zeros(sy);     W3 = zeros(sx,sy);
    dz1 = 0;            dz2 = size(W2W2,2);    dz3 = 0;

    %{
    % Init. shared-latent space.
    [W r] = eig(S, blkdiag(W4W4, W1W1) + Noise);     [r perm] = sort(diag(r), 'descend');
    W = (blkdiag(W4W4, W1W1) + Noise) * W(:,perm(r>1)) * sqrt(diag(r(r>1)-1)); % Retrieve residual variance basis.
    W = W(:, 1:min([sx sy size(W,2)]));   % Upper boundary for dz2 is min(sx,sy).
    W2 = W(sx+1:end,:);     W3 = W(1:sx,:);
    W2W2 = W2*W2';          W3W3 = W3*W3';
    dz2 = size(W,2);
    %}
    disp ' '; diagnostic;

    oldLML = LML_RCA;
    while ~converged
        % RCAs for W1, W2 individually.
        %{
        [W r] = eig(blkdiag(Sx,Sy), blkdiag(W3W3, W2W2) + Noise);   [r perm] = sort(diag(r), 'descend');
        W = (blkdiag(W3W3, W2W2) + Noise) * W(:,perm(r>1)) * sqrt(diag(r(r>1)-1));
        W1 = W(sx+1:end,:);     W4 = W(1:sx,:);
        % W1 = W1(:, 1:min([sy size(W1,2)]));   % Upper boundary for dz1 is sy.
        % W4 = W4(:, 1:min([sx size(W4,2)]));   % Upper boundary for dz3 is sx.
        %}
        [W1 Ly] = eig(Sy, W2W2 + Ny);           [Ly permy] = sort(diag(Ly), 'descend');
        [W4 Lx] = eig(Sx, W3W3 + Nx);           [Lx permx] = sort(diag(Lx), 'descend');
        W1 = (W2W2 + Ny) * W1(:,permy(Ly>1)) * sqrt(diag(Ly(Ly>1)-1));     % Retrieve residual variance basis.
        W4 = (W3W3 + Nx) * W4(:,permx(Lx>1)) * sqrt(diag(Lx(Lx>1)-1));
        W1W1 = W1*W1';                          W4W4 = W4*W4';
        dz1 = size(W1,2);                       dz3 = size(W4,2);
%         disp ' ', disp(['# Iteration ' num2str(k)])
%         disp('## RCA for W1, W2 individually'), diagnostic;

        % RCA for W2, W3 jointly.
        [W r] = eig(S, blkdiag(W4W4, W1W1) + Noise);    [r perm] = sort(diag(r), 'descend');
        W = (blkdiag(W4W4, W1W1) + Noise) * W(:,perm(r>1)) * sqrt(diag(r(r>1)-1)); % Retrieve residual variance basis.
        % Try without boundary!
%         W = W(:, 1:min([sx sy size(W,2)]));   % Upper boundary for dz2 is min(sx,sy). 
        W2 = W(sx+1:end,:);     W3 = W(1:sx,:);
        W2W2 = W2*W2';          W3W3 = W3*W3';
        dz2 = size(W,2);
        disp('## RCA for W2, W3 jointly'), diagnostic;

        %{
        % Assess prediction.
        YTrain_RCA = (W2*W3'*pdinv(W4W4 + W3W3 + Nx)*X')'+ repmat(mean(Y,1), size(X,1), 1);
        Ysqdiff = (Y - YTrain_RCA).^2;       trainingRMSerror_RCA = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff));
        subplot(312), plot(.01*k, RMSerror_RCA,'r.'), % plot(.01*k, trainingRMSerror_RCA,'b.')
        subplot(313), plot(.01*k, LML_RCA,'r.')
        %}
        k = k + 1;
        converged = (LML_RCA - oldLML) < limit;
        oldLML = LML_RCA;
    end
    converged = false;
    stats.alpha(i) = alpha;
    stats.RMSerror(i) = RMSerror_RCA;
    stats.dz1(i) = dz1;
    stats.dz2(i) = dz2;
    stats.dz3(i) = dz3;
    i = i + 1;
end


%% Train with pCCA. minRMS=3.0291 for d=18
for d = 1:57
%     d=18;
    [MLE E var CCA] = pcca(X,Y,d);
    YPred_pCCA = XTest*pdinv(Sx)*MLE.Wzx * MLE.Wzy' + repmat(mean(Y,1), size(XTest,1), 1);
    
    Ysqdiff = (YTest - YPred_pCCA).^2;
    RMSerror_pCCA(d) = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff)); % Root mean-square error.
    % C = [MLE.Wzx; MLE.Wzy]*[MLE.Wzx' MLE.Wzy'] + blkdiag(Sx-(MLE.Wzx*MLE.Wzx'), Sy-(MLE.Wzy*MLE.Wzy'));
    % LML_pCCA = -.5 * ((sy+sx)*log(2*pi) + logdet(C) + trace(sum(sum(pdinv(C)'.*S))));
end
% figure(2), clf, plot(RMSerror_pCCA,'b.'), figure(1);

makeplot


%{
%% Train with CCA.
k = 20; % for k = 1:57
[Wx Wy r U V] = cca(X,Y,k); % Canonical coefficients, correlations & covariates.
Uo = [U ones(n, 1)]; % Bias term for linear regression.
Vo = [V ones(n, 1)];
Wuy = pdinv(Uo'*Uo)*Uo'*Y; % Least squares solution of Uo*Wuy = Y.

UoTest = [XTest*Wx ones(size(XTest,1), 1)];
YPred_CCA = UoTest * Wuy; % + repmat(mean(Y), size(YTest,1), 1);

Ysqdiff = (YTest - YPred_CCA).^2;
RMSerror_CCA = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff)); % Root mean-square error.
% end
% hold on, plot(RMSerror_CCA,'r.');
%}

%% Animate predictions.
demSilhouetteLr1,
%{
clf
xyzankurAnimCompareMultiple(YTest, {YPred_LR, YPred_pCCA, YPred_RCA}, ...
    96, {'Ytest','Linear Regression', ['pCCA (d=' num2str(length(CCA.r)) ')'], 'RCA'});
%}
