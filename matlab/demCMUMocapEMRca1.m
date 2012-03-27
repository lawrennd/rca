% DEMCMUMOCAPEMRCA1 EM-RCA demo on reconstruction of the stick man skeleton
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
importTool({'rca','ndlutil','datasets','mocap'})
addpath(genpath('~/CMUmocap/'))
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
addpath(genpath('~/Desktop/CMUmocap/all_asfamc/subjects/'))
includeMotionCategories

figure(1), clf, colormap('hot');    figure(2), clf, colormap('hot');    figure(3), clf, colormap('hot');
figure(4), clf, colormap('hot');    figure(5), clf, colormap('hot')


Lambda = skelConnectionMatrix(acclaimReadSkel('01.asf'));
d = size(Lambda,1);
Lambda = Lambda + Lambda' + eye(d);
nonZero = find(ones(d));                    % Induce any prior knowledge of zeros. If none, this is all ones.
lambda = 5.^linspace(-8,3,30);
Lambda_hat = cell(length(lambda),1);
GLASSOrocstats = zeros(length(lambda), 4);
EMRCArocstats = zeros(length(lambda), 4); 
options = struct('verbose',1,'order',-1);

figure(1), imagesc(Lambda), colorbar, title('Ground truth network')

%{
[walkSkeletons, walkChannels, walkXYZChannels] = collectSkeletonData(walking, .2, false);
[runSkeletons, runChannels, runXYZChannels] = collectSkeletonData(running, .2, false);
[jumpSkeletons, jumpChannels, jumpXYZChannels] = collectSkeletonData(jumping, .2, false);
[miscSkeletons, miscChannels, miscXYZChannels] = collectSkeletonData({playground, physical_activities_and_sports, situations_and_scenarios }, .2, false);
[allSkeletons, allChannels, allXYZChannels] = collectSkeletonData(categories); % Takes several hours. Run once.
[danceSkeletons, danceChannels, danceXYZChannels] = collectSkeletonData(dance, .1, false);
%}
load walk_jump_run_misc_Channels.mat
walkSkeletons{11} = {}; walkChannels{11} = {}; walkXYZChannels{11} = {}; % remove stealth motions.

%% Compute centred squared-distance matrix (kernel) for every frame.
HKH_sum = zeros(d,d);
H = eye(d) - ones(d)./d;
counter = 0;
Y = zeros(1000000,31);
Yaux = zeros(1000000,31);
rocstats = zeros(1000000,4);
figure(5), handle = xyzVisualise(reshape(jumpXYZChannels{1}{1}(1,:), [], 3), skelConnectionMatrix(jumpSkeletons{1}{1}));

% Channels = {walkChannels runChannels jumpChannels miscChannels}; XYZChannels = {walkXYZChannels runXYZChannels jumpXYZChannels miscXYZChannels};
Channels = {danceChannels}; XYZChannels = {danceXYZChannels};
tic
for iChannel = 1:length(Channels)
    for iSubject = 1:length(Channels{iChannel})
        for iTrial = 1:length(Channels{iChannel}{iSubject})
            for iFrame = 1:ceil(size(Channels{iChannel}{iSubject}{iTrial},1)/1)
                X = reshape(XYZChannels{iChannel}{iSubject}{iTrial}(iFrame,:), [], 3);
%                 X = X + X.*randn(size(X))*1e-1;
                %   xyzModify(handle, X, jumpSkeletons{iSubject}{iTrial});
                HX = H*X;  % HKH = HX*HX';
                Y(counter*3+1:(counter+1)*3,:) = HX';
                Yaux = Y(1:(counter+1)*3,:);	A = Yaux'*Yaux / size(Yaux,1);
                [B, U, jitter] = pdinv(A); % pdinv(HKH_sum);
                if cond(B) < 1e+10
                    counter = counter + 1;
                    %   HKH_sum = HKH_sum + HKH;
                    rocstats(counter,:) = emrcaRocStats(Lambda, B<0);   % Performance evaluation.
                else
                    Y(counter*3+1:(counter+1)*3,:) = zeros(3,d);
                end
            end
        end
    end
end
Y = Y(1:counter*3,:);
Cy = Y'*Y / size(Y,1); % HKH_sum / (counter*3);
jit = randn(counter,1)*.001;  
rocstats = rocstats(1:counter, :);
toc

TPs = rocstats(:,1); FPs = rocstats(:,2); FNs = rocstats(:,3);  Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
figure(2), imagesc(Cy), colorbar, title('Sum of centred sq.distance matrices across frames')
figure(3), clf, hold on, plot(Recalls, Precisions, '.b'), text(Recalls(1:end)+jit, Precisions(1:end)+jit, num2cell(1:counter), 'fontsize',7), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), legend('inv.cov'), title('Recall-Precision')
figure(3), text(Recalls(end), Precisions(end), num2cell(counter), 'fontsize', 8, 'color', 'red', 'fontweight', 'bold')
figure(5), imagesc(pdinv(Cy)), colorbar, title('empirical inverse-covariance ')
%}

%% Standard Glasso on mocap data, with varying lambda.
%
warmLambda_hat =  eye(d); % pdinv(Cy);
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
figure(3), hold on, plot(Recalls, Precisions, '-xr', 'Markersize',10), text(Recalls+jit, Precisions+jit, num2cell(lambda),'fontsize',8), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), legend('Glasso'), title('Recall-Precision')


%% Recovery of sparse-inverse and low-rank covariance via iterative application of EM and RCA.
%
lambda = 5.^linspace(-8,-4,30);
emrca_options = struct('limit',1e-4, 'showProgress',0 , 'verbose',0, 'errorCheck',1, 'maxNumIter',1000);
sigma2_n = .01*trace(Cy);                                   % Noise variance.
[S D] = eig(Cy);     [D perm] = sort(diag(D),'descend');    % Initialise W with a PCA low-rank estimate.
W_hat_old = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));
WWt_hat_old = W_hat_old * W_hat_old';
Lambda_hat_new = cell(length(lambda),1);
Lambda_hat_old = eye(d); % pdinv(Cy);                                 % Initialise Lambda_hat with the empirical inverse-covariance.
parfor (i = 1:length(lambda),8)                             % Try different magnitudes of lambda.
    [WWt_hat_new, Lambda_hat_new{i}, Lambda_hat_new_inv] = emrca(Y, WWt_hat_old, Lambda_hat_old, sigma2_n, lambda(i), nonZero, emrca_options);
    % Plot results.
    figure(5), clf, subplot(121), imagesc(Lambda_hat_new{i}), colorbar, title([ 'EM/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    subplot(122), imagesc(WWt_hat_new), colorbar, title('RCA-recovered WW'''), colorbar
    % Performance stats. Row format in pstats : [ TP FP FN TN ].
    EMRCArocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat_new{i}<0);
end


%% Process performances measures.
TPs = EMRCArocstats(:,1); FPs = EMRCArocstats(:,2); FNs = EMRCArocstats(:,3); TNs = EMRCArocstats(:,4);
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);
%% Plot performance.
figure(3), hold on, plot(Recalls+jit, Precisions+jit, '-o'), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), title('Recall-Precision')
text(Recalls+jit, Precisions+jit, num2cell(lambda)), hold on, legend('inv.cov','Glasso','EM-RCA')
% figure(4), hold on, plot(FPRs, Recalls, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'EM/RCA auc: ' num2str(AUC) ], 4), title('ROC');
%}



