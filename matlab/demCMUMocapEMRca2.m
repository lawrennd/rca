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


clc
addpath(genpath('~/mlprojects/matlab/general/'))
importTool({'rca','ndlutil','mocap'})
addpath(genpath('~/CMUmocap/'))
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
addpath(genpath('~/Desktop/CMUmocap/all_asfamc/subjects/'))
s = RandStream('mcg16807','Seed', 1985); RandStream.setDefaultStream(s) % 1985
includeMotionCategories

figure(1), clf, colormap('hot'); figure(2), clf, colormap('hot'); figure(3), clf, colormap('hot'); figure(4), clf, colormap('hot'); figure(5), clf, colormap('hot')

d = 31;
Lambda = skelConnectionMatrix(acclaimReadSkel('01.asf'));   Lambda = Lambda + Lambda' + eye(d);
nonZero = find(ones(d));                    % Induce any prior knowledge of zeros. If none, this is all ones.
lambda = 5.^linspace(-8,3,30);
Lambda_hat = cell(length(lambda),1);
options = struct('verbose',1,'order',-1);
GLASSOrocstats = zeros(length(lambda), 4);   EMRCArocstats = zeros(length(lambda), 4);
figure(1), imagesc(Lambda), colorbar, title('Ground truth network')

%{
[walkSkeletons, walkChannels, walkXYZChannels] = collectSkeletonData(walking, .01, false);
[danceSkeletons, danceChannels, danceXYZChannels] = collectSkeletonData(dance, .01, false);
%}
load walk_dance_Channels.mat

%% Compute centred squared-distance matrix (kernel) for every frame.
H = eye(d) - ones(d)./d;
counter = 0;
Ytemp = {zeros(100000,31) zeros(100000,31)};
Yaux = zeros(100000,31);
xyz = {zeros(10000,93) zeros(10000,93)};
Channels = {walkChannels danceChannels}; XYZChannels = {walkXYZChannels danceXYZChannels}; Skeletons = {walkSkeletons danceSkeletons};
nFrames = zeros(1,length(Channels));
rocstats = zeros(2,4);
for iCh = 1:length(Channels)
    for iSubject = 1:length(Channels{iCh})
        for iTrial = 1:length(Channels{iCh}{iSubject})
            for iFrame = 1:ceil(size(Channels{iCh}{iSubject}{iTrial},1))
                nFrames(iCh) = nFrames(iCh) + 1;
                xyz{iCh}(nFrames(iCh),:) = XYZChannels{iCh}{iSubject}{iTrial}(iFrame,:);
            end
            %   skelPlayData(Skeletons{iCh}{iSubject}{iTrial}, Channels{iCh}{iSubject}{iTrial})
        end
    end
    isIll = true;
    Ytemp{iCh} = zeros(nFrames(iCh)*3,d);
    while isIll
        nicexyz = xyz{iCh}(randsample(nFrames(iCh),min(3000,nFrames(iCh))),:);
        for i = 1:size(nicexyz,1)
            X = reshape(nicexyz(i,:), [], 3);
            X = X + randn(size(X))*2;
            Ytemp{iCh}((i-1)*3+1:i*3,:) = (H*X)';
        end
        [B, U, jitter] = pdinv(Ytemp{iCh}'*Ytemp{iCh} / size(Ytemp{iCh},1));
        rocstats(iCh,:) = emrcaRocStats(Lambda, B<0);   % Measure of how good B is for initializing Lambda_hat.
        if cond(B) < 1e+10
            isIll = false;
        end
    end
end
%{
% for iCh = 1:length(Channels)
%     for iSubject = 1:length(Channels{iCh})
%         for iTrial = 1:length(Channels{iCh}{iSubject})
%             for iFrame = 1:ceil(size(Channels{iCh}{iSubject}{iTrial},1))
%                 x = XYZChannels{iCh}{iSubject}{iTrial}(iFrame,:);
%                 X = reshape(x, [], 3);
%                 HX = H*X;
%                 Y(nFrames(iCh)*3+1:(nFrames(iCh)+1)*3,:) = HX';
%                 Yaux = Y(1:(nFrames(iCh)+1)*3,:);
%                 A = Yaux'*Yaux / size(Yaux,1);
%                 [B, U, jitter] = pdinv(A);
%                 if cond(B) < 1e+10
%                     nFrames(iCh) = nFrames(iCh) + 1;
%                     xyz{iCh}(nFrames(iCh),:) = XYZChannels{iCh}{iSubject}{iTrial}(iFrame,:);
%                     rocstats(nFrames(iCh),:) = emrcaRocStats(Lambda, B<0);   % Performance evaluation.
%                 else
%                     Y(nFrames(iCh)*3+1:(nFrames(iCh)+1)*3,:) = zeros(3,d);
%                 end
%             end
%             %   skelPlayData(Skeletons{iCh}{iSubject}{iTrial}, Channels{iCh}{iSubject}{iTrial})
%         end
%     end
% end
%}
Y = [Ytemp{1}; Ytemp{2}];   % 1 for walking, 2 for dancing
Cy = Y'*Y / size(Y,1);
figure(2), imagesc(Cy), colorbar, title('Sum of centred sq.distance matrices across frames')
figure(5), imagesc(pdinv(Cy)), colorbar, title('empirical inverse-covariance ')
%}

%% Standard Glasso on mocap data, with varying lambda.
%
warmLambda_hat =  eye(d);%*trace(Cy)/d; % pdinv(Cy);
funObj = @(x)sparsePrecisionObj(x, d, nonZero, Cy);
parfor (i = 1:length(lambda),8)
    Lambda_hat{i} = eye(d);
    Lambda_hat{i}(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda(i)*ones(d*d,1), options);
    % Performance stats. Row format in pstats : [ TP FP FN TN ].
    GLASSOrocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat{i}<0);
    figure(5), imagesc(Lambda_hat{i}), colormap(hot), colorbar, title([ 'GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
end
TPs = GLASSOrocstats(:,1); FPs = GLASSOrocstats(:,2); FNs = GLASSOrocstats(:,3);
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
jit = randn(length(lambda),1)*.001;
figure(3), hold on, plot(Recalls, Precisions, '-xb', 'Markersize',10), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
figure(3), hold on, plot(rocstats(:,1) ./ (rocstats(:,1) + rocstats(:,3)), rocstats(:,1) ./ (rocstats(:,1) + rocstats(:,2)), 'sr'), legend('Glasso','inv.cov')
text(Recalls+jit, Precisions+jit, num2cell(lambda))


%% Recovery of sparse-inverse and low-rank covariance via iterative application of EM and RCA.
lambda = 5.^linspace(-8,-4,30);
emrca_options = struct('limit',1e-4, 'showProgress',0 , 'verbose',0, 'errorCheck',1, 'maxNumIter',250);
sigma2_n = .01*trace(Cy);                                   % Noise variance.
[S D] = eig(Cy);     [D perm] = sort(diag(D),'descend');    % Initialise W with a PCA low-rank estimate.
W_hat_old = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));
WWt_hat_old = W_hat_old * W_hat_old';
Lambda_hat_new = cell(length(lambda),1);
Lambda_hat_old = eye(d); % pdinv(Cy); %                                % Initialise Lambda_hat with the empirical inverse-covariance.
parfor (i = 1:length(lambda),8)                             % Try different magnitudes of lambda.
    [WWt_hat_new, Lambda_hat_new{i}, Lambda_hat_new_inv] = emrca(Y, WWt_hat_old, Lambda_hat_old, sigma2_n, lambda(i), nonZero, emrca_options);
    % Plot results.
    figure(5), clf, subplot(121), imagesc(Lambda_hat_new{i}), colorbar, title([ 'EM/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    subplot(122), imagesc(WWt_hat_new), colorbar, title('RCA-recovered WW'''), colorbar
    % Performance stats. Row format in pstats : [ TP FP FN TN ].
    EMRCArocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat_new{i}<0);
end
%% Plot performance.
TPs = EMRCArocstats(:,1); FPs = EMRCArocstats(:,2); FNs = EMRCArocstats(:,3); TNs = EMRCArocstats(:,4);
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);
figure(3), hold on, plot(Recalls, Precisions, '-og'), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), legend('Glasso','inv.cov','EM-RCA')
text(Recalls+jit, Precisions+jit, num2cell(lambda))
