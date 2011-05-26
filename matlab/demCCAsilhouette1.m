% DEMCCASILHOUETTE Compute canonical variates and correlations on
% silhouette data.
%
% FORMAT
% DESC runs simple CCA on the Agawal and Triggs data.
%
% SEEALSO : cca, pcca
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011

% CCA

addpath(genpath('~/mlprojects/matlab/general/'))
importLatest('datasets');
importLatest('ndlutil');
importLatest('mocap');

clear, clf, clc
dataSetName = 'silhouette';
converged = false;
k=1;

% Load data.
[X, Y, XTest, YTest] = mapLoadData(dataSetName);
[n, sx] = size(X);
sy = size(Y, 2);

Sy = cov(Y, 1); Sx = cov(X,1);  S = cov([X Y], 1);     % Data covariances.
Sx = Sx + abs(mean(diag(Sx)))*1e-6*eye(sx); S(1:sx,1:sx) = Sx;  % Add jitter.


%{
% Interesting initialisations for init-with-CCA variant:
%alpha = .000005   dz2 = 11 : minRMS = 3.0972, maxLML = -86.296
%alpha = .0000005  dz2 = 11 : minRMS = 3.0361, maxLML = -73.6919
%alpha = .0000002  dz2 = 11 : minRMS = 3.0577, maxLML = -69.3348
%alpha = .0000002  dz2 = 13 : minRMS = 3.0295, maxLML = -68.1126
%alpha = .0000002  dz2 = 14 : minRMS = 3.0473, maxLML = -67.5159
%alpha = .0000002  dz2 = 15 : minRMS = 3.0729, maxLML = -66.8512
%alpha = .0000002  dz2 = 16 : minRMS = 3.0673, maxLML = -66.3111
%alpha = .0000002  dz2 = 17 : minRMS = 3.0383, maxLML = -65.7893
%alpha = .0000002  dz2 = 18 : minRMS = 3.0375, maxLML = -65.6238
%alpha = .0000002  dz2 = 19 : minRMS = 3.0614, maxLML = -65.1887
%alpha = .0000002  dz2 = 20 : minRMS = 3.0071, maxLML = -64.3489
%alpha = .0000002  dz2 = 21 : minRMS = 3.0096, maxLML = -64.263
%alpha = .0000002  dz2 = 22 : minRMS = 3.0089, maxLML = -63.8763
%alpha = .0000002  dz2 = 23 : minRMS = 3.032,  maxLML = -62.9842
%alpha = .0000002  dz2 = 24 : minRMS = 3.0342, maxLML = -62.9466
%alpha = .0000002  dz2 = 25 : minRMS = 3.0692, maxLML = -62.5618
%alpha = .0000002  dz2 = 26 : minRMS = 3.0913, maxLML = -62.2339
%alpha = .0000002  dz2 = 27 : minRMS = 3.0913, maxLML = -62.2339

%}

% Fit private and shared latent spaces with RCA-CCA.
% Initialise noise variances.
alpha = .0000002;  % Fraction of residual variance (removed pCCA-shared variance) as noise.
dz2 = 11;
limit = 1e-5;

%{
% Init. private-latent spaces with PPCA.
[~, Ly] = eig(Sy);  [Ly, ~] = sort(diag(Ly),'descend');
[~, Lx] = eig(Sx);  [Lx, ~] = sort(diag(Lx),'descend');
ny = .05*mean(Ly);    nx = .05*mean(Lx);
Ny = ny*eye(sy);        Nx = nx*eye(sx);        Noise = blkdiag(Nx, Ny);
% Ly = Ly(Ly > ny);       Lx = Lx(Lx > nx);
Ly = Ly(cumsum(Ly) < (1-alpha)*trace(Sy));      Lx = Lx(cumsum(Lx) < (1-alpha)*trace(Sx));
Ly = Ly(cumsum(Ly) <= (1-beta)*sum(Ly));        Lx = Lx(cumsum(Lx) <= (1-beta)*sum(Lx));
dz1 = length(Ly);       dz3 = length(Lx);
Uy = Uy(:, Iy(1:dz1));  Ux = Ux(:, Ix(1:dz3));
W1 = Uy*sqrt(diag(Ly)-ny*eye(dz1)); % ML estimates of PPCA latent bases, given desired noise percentage.
W4 = Ux*sqrt(diag(Lx)-nx*eye(dz3));
W1W1 = W1*W1';          W4W4 = W4*W4';
%}

%%{
% Init. shared-latent space with PCCA.
MLE = pcca(X,Y,dz2);
W2 = MLE.Wzy;                   W3 = MLE.Wzx;      % ML estimates of pCCA mappings.
W2W2 = W2*W2';                  W3W3 = W3*W3';
oldtraceW2W2 = trace(W2W2);     oldtraceW3W3 = trace(W3W3);  % Track traces of W2W2, W3W3.
PSIy_aux = Sy - W2W2;           PSIx_aux = Sx - W3W3;   % Residual covariance.
[~, Ly] = eig(PSIy_aux);        [Ly, ~] = sort(diag(Ly),'descend'); % Eigenvalue analysis on residual covariances.
[~, Lx] = eig(PSIx_aux);        [Lx, ~] = sort(diag(Lx),'descend');
Ly = max(Ly, 0);                Lx = max(Lx, 0);    % Remove round-off errors.
% ny = alpha*trace(PSIy_aux)/sy;  nx = alpha*trace(PSIx_aux)/sx;
posy = cumsum(Ly)>(1-alpha)*trace(PSIy_aux);    posx = cumsum(Lx)>(1-alpha)*trace(PSIx_aux);
ny = sum(Ly(posy(1:end)))/sy; nx = sum(Lx(posx(1:end)))/sx;
% dz1 = sy;                       dz3 = sx;
dz1 = sum(~posy);               dz3 = sum(~posx);
Ny = ny*eye(sy);                Nx = nx*eye(sx);        Noise = blkdiag(Nx, Ny);
PSIy = PSIy_aux - Ny;           PSIx = PSIx_aux - Nx;   % Also remove noises from residuals.
 
        % Diagnostic: Checking the RMS error and LML of the model so far.
            % YPred_pCCA = XTest * CCA.Wx * sqrt(diag(CCA.r)) * MLE.Wzy' + repmat(mean(Y,1), size(XTest,1), 1); % graphical way
        YPred_pCCA = (W2*W3'*pdinv(Sx)*XTest')'+ repmat(mean(Y,1), size(XTest,1), 1); % conditional-Gaussian way
        Ysqdiff = (YTest - YPred_pCCA).^2;
        RMSerror_pCCA = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff));
        C = [W3;W2]*[W3' W2'] + Noise;      % restriced pCCA covariance model (thrown away parts not explained by W2,W3,noise)
        LML_pCCA = -.5 * ((sy+sx)*log(2*pi) + logdet(C) + trace(sum(sum(pdinv(C)'.*S))));
        disp ' ', disp('## pCCA initialisation')
        disp(['####    Check trace balance: traceC = ' num2str(trace(C)) '    traceS = ' num2str(trace(S))])
        disp(['####    LML_pCCA = ' num2str(LML_pCCA)])
        disp(['####    RMSerror_pCCA = ' num2str(RMSerror_pCCA)]);

oldLML = LML_pCCA;
while ~converged
    % Two RCAs for W1, W2 individually.
    [W1, Ly] = eig(Sy, W2W2 + Ny);     [Ly, Iy] = sort(diag(Ly), 'descend');
    [W4, Lx] = eig(Sx, W3W3 + Nx);     [Lx, Ix] = sort(diag(Lx), 'descend');
    Ly = max(Ly, 0);                        Lx = max(Lx, 0);    % Remove round-off errors.
    % Retrieve residual basis.
    W1 = W1(:, Iy(1:dz1));                  W4 = W4(:, Ix(1:dz3));
    W1 = (W2W2 + Ny) * W1 * sqrt(diag(Ly(1:dz1))-eye(dz1));
    W4 = (W3W3 + Nx) * W4 * sqrt(diag(Lx(1:dz3))-eye(dz3));
    W1W1 = W1*W1';                          W4W4 = W4*W4';
%     W1W1 = W1W1*trace(PSIy)/trace(W1W1);    W4W4 = W4W4*trace(PSIx)/trace(W4W4);
    
        %%{
        % Diagnostic: Checking the RMS error and LML of the model so far.
        C = blkdiag(W4W4,W1W1) + [W3;W2]*[W3' W2'] + Noise;
        LML_CcaRca = -.5 * ((sy+sx)*log(2*pi) + logdet(C) + trace(sum(sum(pdinv(C)'.*S))));
%         subplot(221), imshow(C), subplot(222), imshow(S)
        YPred_CcaRca = (W2*W3'*pdinv(W4W4 + W3W3 + Nx)*XTest')'+ repmat(mean(Y,1), size(XTest,1), 1);
        Ysqdiff = (YTest - YPred_CcaRca).^2; RMSerror_CcaRca = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff));
        disp ' ', disp(['# Iteration ' num2str(k)]), disp('## RCA for W1, W2 individually')
        disp(['#### Check trace balance: traceC = ' num2str(trace(C)) '    traceS = ' num2str(trace(S))])
        disp(['#### LML_CcaRca = ' num2str(LML_CcaRca)])
        disp(['#### RMSerror_CcaRca = ' num2str(RMSerror_CcaRca)])
        %}    
    
    % RCA for W2, W3 jointly.
    [W r] = eig(S, blkdiag(W4W4, W1W1) + Noise);            % RCA form.
    [r, I] = sort(diag(r), 'descend');  r = max(r, 0);      % Remove round-off errors.
    % Retrieve residual basis.
    Wy = W(sx+1:end, I(1:dz2));             Wx = W(1:sx, I(1:dz2));
    W2 = (W1W1 + Ny) * Wy * sqrt(diag(r(1:dz2))-eye(dz2));
    W3 = (W4W4 + Nx) * Wx * sqrt(diag(r(1:dz2))-eye(dz2));
%     if k==1
%         W2 = W2*sqrt(oldtraceW2W2/trace(W2*W2')); W3 = W3*sqrt(oldtraceW3W3/trace(W3*W3'));
%     end
    W2W2 = W2*W2';                          W3W3 = W3*W3';
%     oldtraceW2W2 = trace(W2W2);             oldtraceW3W3 = trace(W3W3);  % Update traces of W2W2, W3W3.

    %{
    % Two RCAs for W1, W2 individually.
    [W1, Ly] = eig(Sy, W2W2 + Ny);      [Ly, Iy] = sort(diag(Ly), 'descend');
    [W4, Lx] = eig(Sx, W3W3 + Nx);      [Lx, Ix] = sort(diag(Lx), 'descend');
    Ly = max(Ly, 0);                    Lx = max(Lx, 0);    % Remove round-off errors.
    % Retrieve residual basis.
    W1 = W1(:, Iy(1:dz1));              W4 = W4(:, Ix(1:dz3));
    W1 = (W2W2 + Ny) * W1 * sqrt(diag(Ly(1:dz1)));
    W4 = (W3W3 + Nx) * W4 * sqrt(diag(Lx(1:dz3)));
    W1W1 = W1*W1';                      W4W4 = W4*W4';
    %}

        %%{
        % Diagnostic: Checking the RMS error and LML of the model so far.
        C = (blkdiag(W4W4,W1W1) + [W3;W2]*[W3' W2'] + Noise);
        LML_CcaRca = -.5 * ((sy+sx)*log(2*pi) + logdet(C) + trace(sum(sum(pdinv(C).*S))));
%         subplot(221), imshow(C), subplot(222), imshow(S)
        YPred_CcaRca = (W2*W3'*pdinv(W4W4 + W3W3 + Nx)*XTest')'+ repmat(mean(Y,1), size(XTest,1), 1);
        Ysqdiff = (YTest - YPred_CcaRca).^2; RMSerror_CcaRca = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff));
        disp('## RCA for W2, W3 jointly')
        disp(['#### Check trace balance: traceC = ' num2str(trace(C)) '    traceS = ' num2str(trace(S))])
        disp(['#### LML_CcaRca = ' num2str(LML_CcaRca)])
        disp(['#### RMSerror_CcaRca = ' num2str(RMSerror_CcaRca)])
        %}

        % Assess prediction.
        YTrain_CcaRca = (W2*W3'*pdinv(W4W4 + W3W3 + Nx)*X')'+ repmat(mean(Y,1), size(X,1), 1);
        Ysqdiff = (Y - YTrain_CcaRca).^2;       trainingRMSerror_CcaRca = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff));
        subplot(121), hold on, plot(.01*k, RMSerror_CcaRca,'r.'), plot(.01*k, trainingRMSerror_CcaRca,'b.'), %title('RMS error');
        subplot(122), hold on, plot(.01*k, LML_CcaRca,'r.') %title('Log-marginal likelihood');
        k = k + 1;
        
    converged = (LML_CcaRca - oldLML) < limit;
    oldLML = LML_CcaRca;
end

    
%% Train with pCCA.
d = dz2; % for k = 1:57
[MLE E var CCA] = pcca(X,Y,d);
ZTest = XTest * CCA.Wx * sqrt(diag(CCA.r)); % E(Zt|Xt).

YPred_pCCA = ZTest * MLE.Wzy' + repmat(mean(Y,1), size(ZTest,1), 1);
Ysqdiff = (YTest - YPred_pCCA).^2;
RMSerror_pCCA = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff)); % Root mean-square error.

C = [MLE.Wzx; MLE.Wzy]*[MLE.Wzx' MLE.Wzy'] + blkdiag(Sx-(MLE.Wzx*MLE.Wzx'), Sy-(MLE.Wzy*MLE.Wzy'));
LML_pCCA = -.5 * ((sy+sx)*log(2*pi) + logdet(C) + trace(sum(sum(pdinv(C)'.*S))));
% end
% plot(RMSerror_pCCA,'b.');


%% Train with CCA.
k = 57; % for k = 1:57
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


%% Animate predictions.
demLinRegSilhouette, %clf
% xyzankurAnimCompare(YTest, YPred_LR, 96, {'Ytest','Linear Regression',...
%     ['CCA (d=' num2str(k) ')'], ['pCCA (d=' num2str(length(CCA.r)) ')'], 'CcaRca'},...
%     YPred_CCA, YPred_pCCA, YPred_CcaRca);
