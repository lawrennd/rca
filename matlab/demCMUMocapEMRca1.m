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
addpath(genpath('~/mlprojects/rca/matlab/glasso/'))
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
importTool({'rca','ndlutil','datasets','mocap'})
addpath(genpath('~/Desktop/CMUmocap/all_asfamc/subjects/'))
includeMotionCategories

figure(1), clf, colormap('hot');    figure(2), clf, colormap('hot')
figure(3), clf, colormap('hot');    figure(4), clf, colormap('hot')
figure(5), clf, colormap('hot')


d = 31;
Lambda = skelConnectionMatrix(acclaimReadSkel('01.asf'));   Lambda = Lambda + Lambda' + eye(d);
nonZero = find(ones(d));                    % Induce any prior knowledge of zeros. If none, this is all ones.
Lambda_hat = cell(length(lambda),1);
rocstats = zeros(length(lambda), 4);
Cor = zeros(2068, 1);
showProgress = 0;

figure(1), imagesc(Lambda), colorbar

%{
[walkSkeletons, walkChannels, walkXYZChannels] = collectSkeletonData(walking, .2, false);
[runSkeletons, runChannels, runXYZChannels] = collectSkeletonData(running, .2, false);
[jumpSkeletons, jumpChannels, jumpXYZChannels] = collectSkeletonData(jumping, .2, false);
[miscSkeletons, miscChannels, miscXYZChannels] = collectSkeletonData({playground, physical_activities_and_sports, situations_and_scenarios }, .2, false);
[allSkeletons, allChannels, allXYZChannels] = collectSkeletonData(categories); % Takes several hours. Run once.
%}

%
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
            D(counter,:) = reshape(ones(d,1)*diag(X*X')' - 2*(X*X') + diag(X*X')*ones(1,d) , 1,[]);
            [B, U, jitter] = pdinv(HKH_sum);
%             figure(5), imagesc(abs(B)); title(num2str(counter)), colorbar
            rocstats(counter,:) = emrcaRocStats(Lambda, (B<0), 0);   % Performance evaluation.
        end
    end
end
DistVar = reshape(var(D), d,[]);
toc

TPs = rocstats(:,1); FPs = rocstats(:,2); FNs = rocstats(:,3); TNs = rocstats(:,4); 
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);
figure(3), clf, hold on, plot(Recalls, Precisions, '.b'), text(Recalls, Precisions, num2cell(1:counter), 'fontsize',7), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), legend('raw'), title('Recall-Precision')
figure(4), clf, hold on, plot(FPRs, Recalls, '.b'), xlim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'GL auc: ' num2str(AUC) ], 4), title('ROC');

% HKH_sum = HKH_sum + eye(d)*trace(HKH_sum)*1e-6;
figure(2), imagesc(HKH_sum), colorbar
figure(5), imagesc(pdinv(HKH_sum)), colorbar
finalSum_ROCstats = emrcaRocStats(Lambda, pdinv(HKH_sum)<0, 0);   % Performance evaluation.
TP = finalSum_ROCstats(1); FP = finalSum_ROCstats(2); FN = finalSum_ROCstats(3); TN = finalSum_ROCstats(4);
Recall = TP / (TP + FN);   Precision = TP / (TP + FP);
finalSum_ROCstats = emrcaRocStats(Lambda, pdinv(HKH_sum)<0, 0);   % Performance evaluation.
TP = finalSum_ROCstats(1); FP = finalSum_ROCstats(2); FN = finalSum_ROCstats(3); TN = finalSum_ROCstats(4);
Recall = TP / (TP + FN);   Precision = TP / (TP + FP);

figure(3), plot(Recall, Precision, 'xr', 'MarkerSize', 10, 'LineWidth',2)
%}


%% Standard Glasso on mocap data, with varying lambda.
lambda = 5.^linspace(5,15,30);
rocstats = zeros(length(lambda), 4);
warmLambda_hat = pdinv(HKH_sum); % eye(d);
funObj = @(x)sparsePrecisionObj(x, d, nonZero, HKH_sum);
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

figure(3), hold on, plot(Recalls, Precisions, '-xr'), text(Recalls, Precisions, num2cell(lambda),'fontsize',8), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), legend('Glasso'), title('Recall-Precision')
figure(4), hold on, plot(FPRs, Recalls, '-xb'), xlim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'GL auc: ' num2str(AUC) ], 4), title('ROC');
%}


