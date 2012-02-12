% DEMCMUMOCAPEMRCA1 EM-RCA demo on reconstruction of the stick man skeleton
% based on 3-D sensor data from motions across the CMU mocap database.
%
% FORMAT
% DESC
%
% SEEALSO : includeMotionCategories
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011, 2012
%
% RCA


clc
addpath(genpath('~/mlprojects/matlab/general/'))
importTool({'rca','ndlutil','datasets','mocap'})
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
addpath(genpath('~/Desktop/CMUmocap/all_asfamc/subjects/'))
includeMotionCategories

figure(1), clf, colormap('hot');    figure(2), clf, colormap('hot')
figure(3), clf, colormap('hot');    figure(4), clf, colormap('hot')
figure(5), clf, colormap('hot')


d = 31;
Lambda = skelConnectionMatrix(acclaimReadSkel('01.asf'));   Lambda = Lambda + Lambda' + eye(d);
nonZero = find(ones(d));                    % Induce any prior knowledge of zeros. If none, this is all ones.
% lambda = 5.^linspace(5,15,30);
lambda = 5.^linspace(1,5,30);
Lambda_hat = cell(length(lambda),1);
rocstats = zeros(length(lambda), 4);
Cor = zeros(2068, 1);
showProgress = 0;

figure(1), imagesc(Lambda), colorbar, title('Ground truth network')

%{
[walkSkeletons, walkChannels, walkXYZChannels] = collectSkeletonData(walking, .2, false);
[runSkeletons, runChannels, runXYZChannels] = collectSkeletonData(running, .2, false);
[jumpSkeletons, jumpChannels, jumpXYZChannels] = collectSkeletonData(jumping, .2, false);
[miscSkeletons, miscChannels, miscXYZChannels] = collectSkeletonData({playground, physical_activities_and_sports, situations_and_scenarios }, .2, false);
[allSkeletons, allChannels, allXYZChannels] = collectSkeletonData(categories); % Takes several hours. Run once.
%}

% Compute centred squared-distance matrix for every frame.
HKH_sum = zeros(d,d);
H = eye(d) - ones(d)./d;
counter = 0;
tic
% handle = xyzVisualise(reshape(jumpXYZChannels{1}{1}(1,:), [], 3), jumpSkeletons{1}{1});

for i = 1:length(jumpChannels)
    for j = 1:length(jumpChannels{i})
        for iFrame = 1:ceil(size(jumpChannels{i}{j},1)/10)
            X = reshape(jumpXYZChannels{i}{j}(iFrame,:), [], 3);
%             pause(1/120); xyzModify(handle, X, jumpSkeletons{i}{j});
            HKH = H*(X*X')*H;
            HKH_sum = HKH_sum + HKH;
            counter = counter + 1;
%             D(counter,:) = reshape(ones(d,1)*diag(X*X')' - 2*(X*X') + diag(X*X')*ones(1,d) , 1,[]);
            [B, U, jitter] = pdinv(HKH_sum);
%             figure(5), imagesc(abs(B)); title(num2str(counter)), colorbar
            rocstats(counter,:) = emrcaRocStats(Lambda, (B<0), 0);   % Performance evaluation.
        end
    end
end
% DistVar = reshape(var(D), d,[]); figure(2), imagesc(DistVar), colorbar
toc

TPs = rocstats(:,1); FPs = rocstats(:,2); FNs = rocstats(:,3);  Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
figure(3), jit = randn(counter,1)*.001; clf, hold on, plot(Recalls, Precisions, '.b'), text(Recalls(1:end)+jit, Precisions(1:end)+jit, num2cell(1:counter), 'fontsize',7), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), legend('inv.cov'), title('Recall-Precision')
% TNs = rocstats(:,4); FPRs = FPs ./ (FPs + TNs);    AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);   figure(4), clf, hold on, plot(FPRs, Recalls, '.b'), xlim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'GL auc: ' num2str(AUC) ], 4), title('ROC');

C = HKH_sum ./ counter;
figure(2), imagesc(C), colorbar, title('Sum of centred sq.distance matrices across frames')
figure(5), imagesc(pdinv(C)), colorbar, title('empirical inverse-covariance ')
% roc = emrcaRocStats(Lambda, pdinv(C)<0, 0);   % Performance evaluation.
% figure(3), plot(roc(1)/(roc(1)+roc(3)), roc(1)/(roc(1)+roc(2)), 'xr', 'MarkerSize', 10, 'LineWidth',2)
% figure(3), plot(rocstats(end,1)/(rocstats(end,1)+rocstats(end,3)), rocstats(end,1)/(rocstats(end,1)+rocstats(end,2)), '.r', 'MarkerSize', 10, 'LineWidth',2)
figure(3), text(Recalls(end), Precisions(end), num2cell(counter), 'fontsize', 8, 'color', 'red', 'fontweight', 'bold')
%}


%% Standard Glasso on mocap data, with varying lambda.
rocstats = zeros(length(lambda), 4);
warmLambda_hat =  pdinv(C);
funObj = @(x)sparsePrecisionObj(x, d, nonZero, C);
options.verbose = 1;    options.order = -1;
for i = 1:length(lambda)
    Lambda_hat{i} = eye(d);
    Lambda_hat{i}(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda(i)*ones(d*d,1), options);
    rocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat{i}<0, 0);   % Performance evaluation.
    
    figure(5), imagesc(Lambda_hat{i}), colormap(hot), colorbar, title([ 'GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
end
TPs = rocstats(:,1); FPs = rocstats(:,2); FNs = rocstats(:,3); TNs = rocstats(:,4); 
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);

jit = randn(length(lambda),1)*.001;
figure(3), hold on, plot(Recalls, Precisions, '-xr'), text(Recalls+jit, Precisions+jit, num2cell(lambda),'fontsize',8), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), legend('Glasso'), title('Recall-Precision')
% figure(4), hold on, plot(FPRs, Recalls, '-xb'), xlim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'GL auc: ' num2str(AUC) ], 4), title('ROC');



