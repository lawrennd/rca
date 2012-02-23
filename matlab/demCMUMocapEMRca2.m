% DEMCMUMOCAPEMRCA2 EM-RCA demo on reconstruction of the stick man skeleton
% based on 3-D sensor data from motions across the CMU mocap database.
%
% FORMAT
% DESC
%
% SEEALSO : includeMotionCategories, emrca
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011, 2012
%
% RCA

% clc
addpath(genpath('~/mlprojects/matlab/general/'))
importTool({'rca','ndlutil','mocap'})
addpath(genpath('~/CMUmocap/'))
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
addpath(genpath('~/Desktop/CMUmocap/all_asfamc/subjects/'))
% s = RandStream('mcg16807','Seed', 1e+6); RandStream.setDefaultStream(s)
includeMotionCategories

% figure(1), clf, colormap('hot'); figure(2), clf, colormap('hot'); figure(5), clf, colormap('hot')

% General settings.
Lambda = skelConnectionMatrix(acclaimReadSkel('01.asf'));
d = size(Lambda,1);
Lambda = Lambda + Lambda' + eye(d);
nonZero = find(ones(d));                    % Induce any prior knowledge of zeros. If none, this is all ones.
lambda = 5.^linspace(-8,3,30);
options = struct('verbose',0,'order',-1);
GLASSOrocstats = zeros(length(lambda), 4);
EMRCArocstats = zeros(length(lambda), 4);
% figure(1), imagesc(Lambda), colorbar, title('Ground truth network')

%{
[walkSkeletons, walkChannels, walkXYZChannels] = collectSkeletonData(walking, .01, false);
[jumpSkeletons, jumpChannels, jumpXYZChannels] = collectSkeletonData(jumping, .1, false);
[runSkeletons, runChannels, runXYZChannels] = collectSkeletonData(running, .1, false);
[danceSkeletons, danceChannels, danceXYZChannels] = collectSkeletonData(dance, .1, false);
%}
load walk_jump_run_dance_UniformlySampled_Channels.mat
walkSkeletons{11} = {}; walkChannels{11} = {};  walkXYZChannels{11} = {};           % Remove stealth motions.
jumpSkeletons{3} = {};  jumpChannels{3} = {};   jumpXYZChannels{3} = {};            % Remove hoping on one foot.

%% Compute centred squared-distance matrix (kernel) for every frame.
H = eye(d) - ones(d)./d;
rocstats = zeros(1,4);
Channels = {walkChannels jumpChannels runChannels danceChannels};
Ytemp = cell(length(Channels), 1);
%{
Ytemp = {zeros(100000,31) zeros(100000,31)};
Yaux = zeros(100000,31);
xyz = {zeros(10000,93) zeros(10000,93)};
XYZChannels = {walkXYZChannels jumpXYZChannels runXYZChannels danceXYZChannels};
Skeletons = {walkSkeletons jumpSkeletons runSkeletons danceSkeletons};
nFrames = zeros(1,length(Channels));
for iCh = 1:length(Channels)
    for iSubject = 1:length(Channels{iCh})
        for iTrial = 1:length(Channels{iCh}{iSubject})
            for iFrame = 1:ceil(size(Channels{iCh}{iSubject}{iTrial},1))
                nFrames(iCh) = nFrames(iCh) + 1;
                xyz{iCh}(nFrames(iCh),:) = XYZChannels{iCh}{iSubject}{iTrial}(iFrame,:);
            end
            % if (iCh > 1)
            %   skelPlayData(Skeletons{iCh}{iSubject}{iTrial}, Channels{iCh}{iSubject}{iTrial})
            % end
        end
    end
    isIll = true;
    Ytemp{iCh} = zeros(nFrames(iCh)*3,d);
    while isIll
        nicexyz = xyz{iCh}(randsample(nFrames(iCh),min(2000,nFrames(iCh))),:);
        for i = 1:size(nicexyz,1)
            X = reshape(nicexyz(i,:), [], 3);
            X = X + randn(size(X))*1;
            Ytemp{iCh}((i-1)*3+1:i*3,:) = (H*X)';
        end
        [B, U, jitter] = pdinv(Ytemp{iCh}'*Ytemp{iCh} / size(Ytemp{iCh},1));
        rocstats(iCh,:) = emrcaRocStats(Lambda, B<0);   % Measure of how good B is for initializing Lambda_hat.
        if cond(B) < 1e+10
            isIll = false;
        end
    end
end
%}
isIll = true;
while isIll
    for iCh = 1:length(Channels)
        Ytemp{iCh} = zeros(nFrames(iCh)*3,d);
        nicexyz = xyz{iCh}(randsample(nFrames(iCh), ceil(nFrames(iCh)*.9)), :);     % Sub-sample from data for stability selection.
        for i = 1:size(nicexyz,1)
            X = reshape(nicexyz(i,:), [], 3);
            X = X + randn(size(X))*1;               % Induce observation noise.
            Ytemp{iCh}((i-1)*3+1:i*3,:) = (H*X)';
        end
    end
    Y = [Ytemp{1}; Ytemp{2}; Ytemp{3}; Ytemp{4}];                                   % 1 for walking, 2 for jumping, 3 for running, 4 for dancing.
    [B, U, jitter] = pdinv(Y'*Y / size(Y,1));
    rocstats = emrcaRocStats(Lambda, B<0);                                          % Measure of how good B is for initializing Lambda_hat.
    if cond(B) < 1e+10
        isIll = false;
    end
end
% Y = Ytemp{2};   % 1 for walking, 2 for dancing
Cy = Y'*Y / size(Y,1);
% figure(2), imagesc(Cy), colorbar, title('Sum of centred sq.distance matrices across frames')
% figure(5), imagesc(pdinv(Cy)), colorbar, title('empirical inverse-covariance ')


%% Standard Glasso on mocap data, with varying lambda.
%
Lambda_hat_glasso = cell(length(lambda),1);
warmLambda_hat =  eye(d); 
funObj = @(x)sparsePrecisionObj(x, d, nonZero, Cy);
parfor (i = 1:length(lambda),8)
    Lambda_hat_glasso{i} = eye(d);
    Lambda_hat_glasso{i}(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda(i)*ones(d*d,1), options);
    % Performance stats. Row format in pstats : [ TP FP FN TN ].
    GLASSOrocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat_glasso{i}<0);
%     figure(5), imagesc(Lambda_hat_glasso{i}), colormap(hot), colorbar, title([ 'GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
end
%{
TPs = GLASSOrocstats(:,1); FPs = GLASSOrocstats(:,2); FNs = GLASSOrocstats(:,3);
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
jit = randn(length(lambda),1)*.001;
figure(2), hold on, plot(Recalls, Precisions, '-xb', 'Markersize',10), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
figure(2), plot(rocstats(:,1) ./ (rocstats(:,1) + rocstats(:,3)), rocstats(:,1) ./ (rocstats(:,1) + rocstats(:,2)), 'sr'), legend('Glasso','inv.cov')
% text(Recalls+jit, Precisions+jit, num2cell(lambda))
%}

%% Recovery of sparse-inverse and low-rank covariance via iterative application of EM and RCA.
% lambda = 5.^linspace(-8,-4,30);
emrca_options = struct('limit',1e-4, 'showProgress',0 , 'verbose',0, 'errorCheck',1, 'maxNumIter',1000);
sigma2_n = 0.5*trace(Cy)/d; % trace(Cy)*1e-3;                                   % Noise variance.
[S D] = eig(Cy);    [D perm] = sort(diag(D),'descend');     % Initialise W with a PCA low-rank estimate.
W_hat_old = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));
WWt_hat_old = W_hat_old * W_hat_old';
Lambda_hat_emrca = cell(length(lambda),1);
Lambda_hat_old = eye(d);
parfor (i = 1:length(lambda),8)                             % Try different magnitudes of lambda.
    [WWt_hat_new, Lambda_hat_emrca{i}] = emrca(Y, WWt_hat_old, Lambda_hat_old, sigma2_n, lambda(i), nonZero, emrca_options);
    % Plot results.
%     figure(5), clf, subplot(121), imagesc(Lambda_hat_emrca{i}), colorbar, title([ 'EM/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
%     subplot(122), imagesc(WWt_hat_new), colorbar, title('RCA-recovered WW'''), colorbar
    % Performance stats. Row format in pstats : [ TP FP FN TN ].
    EMRCArocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat_emrca{i}<0);
end
%% Plot performance.
%{
TPs = EMRCArocstats(:,1); FPs = EMRCArocstats(:,2); FNs = EMRCArocstats(:,3); TNs = EMRCArocstats(:,4);
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);
figure(2), hold on, plot(Recalls, Precisions, '-og'), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), legend('Glasso','inv.cov','EM-RCA')
% text(Recalls+jit, Precisions+jit, num2cell(lambda))
%}