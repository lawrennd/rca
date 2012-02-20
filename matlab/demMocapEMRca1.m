% DEMMOCAPEMRCA1 EM-RCA demo on reconstruction of the stick man
% from a single motion capture file.
%
% FORMAT
% DESC
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011, 2012
%
% RCA

clc
addpath(genpath('~/mlprojects/matlab/general/'))
importTool({'rca','ndlutil','datasets','mocap'})
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
s = RandStream('mcg16807','Seed', 1985); RandStream.setDefaultStream(s) % 1985

figure(1), clf, colormap('hot');    figure(2), clf, colormap('hot')
figure(3), clf, colormap('hot');    figure(4), clf, colormap('hot')
figure(5), clf, colormap('hot')

lambda = 5.^linspace(-8,-5,10);
Lambda_hat = cell(length(lambda),1);    % Sigma_hat = cell(length(lambda),1); 
GLASSOrocstats = zeros(length(lambda), 4);   EMRCArocstats = zeros(length(lambda), 4);
emrca_options = struct('limit',1e-4, 'showProgress',0 , 'verbose',0, 'errorCheck',0, 'maxNumIter',1000);
options = struct('verbose',1,'order',-1);

%% Import mocap data.
[points, pointNames] = mocapParseText('walkSitJog.txt');
connect = mocapConnections('connections_walkJogRun.txt', pointNames);
Y0 = mocapLoadTextData('walkSitJog', false);    Y = Y0(:,1:62); % or stick
[n,d] = size(Y);
connect(5,4) = 1; connect(16,15) = 1; connect(17,16) = 1;   % Add a few missing (reasonable) edges.
Lambda = kron(eye(2), connect + connect' + eye(31));       % Ground truth Lambda.
    %{
        [~, Y0, ~, Ytest0] = mapLoadData('silhouette');
        Y = Y0(:,4:end); Ytest = Ytest0(:,4:end); % Skip the first sensor.
        Y(:,[1:3:end 2:3:end 3:3:end]) = [Y(:,3:3:end) -Y(:,1:3:end) -Y(:,2:3:end)]; % Every 3 columns represent x-y-z coordinates.
        Ytest(:,[1:3:end 2:3:end 3:3:end]) = [Ytest(:,3:3:end) -Ytest(:,1:3:end) -Ytest(:,2:3:end)];
        joint = reshape(Y(1,:)',3,[])';
        % Organise stick man limbs.
        limb{1} = [1 2; 2 3; 3 4];                  % Spine
        limb{2} = [2 5; 5 6; 6 7; 7 8];             % Left-arm
        limb{3} = [2 9; 9 10; 10 11; 11 12];        % Right-arm
        limb{4} = [1 13; 13 14; 14 15; 15 19];      % Left-leg
        limb{5} = [1 16; 16 17; 17 18; 18 20];      % Right-leg
        alllimbs = [limb{1}; limb{2}; limb{3}; limb{4}; limb{5}];
        % Lambda = full( sparse([alllimbs(:,1)', 1:20], [alllimbs(:,2)',1:20], ones(1,39)));
        Yhete = Y;                                  % Training set contains heterogenous movement (confounded dataset).
        Yhomo = Ytest;                              % Test set contains homogenous movement.
        % *** From here on, Y is denoted as the dataset of choice. ***
        Y = Ytest(:,3:3:end);                       % Use z coordinates.
    %}
% selection = randperm(size(Y,1));
% Y = Y(selection(1:fix(size(Y,1)/1)),:);     % Use 100% of the dataset.
nonZero = find(ones(d));                    % Induce any prior knowledge of zeros. If none, this is all ones.
% Normalisation.
Y = Y - repmat(mean(Y),n,1);
stdY = std(Y); stdY(stdY == 0) = 1;
Y = Y ./ repmat(stdY,n,1);
Cy = Y' * Y / n;                            % Covariance estimation.

% Visualisation.
figure(1),subplot(121),
    % axis equal, set(gca,'XLim',[-15 15],'YLim',[-15 15],'ZLim',[0 70]), view(3), grid on,
    % xyzankurDraw([joint(1,:); joint]); text(joint(:,1), joint(:,2), joint(:,3), num2cell(1:20), 'fontsize', 14, 'fontweight', 'bold', 'color', 'r')
handle = stickVisualise(Y0(1, :), connect);  xlim([-.05 .1]), ylim([-.05 .05]), zlim([0 .1])
text(Y0(1,1:31),Y0(1,32:62),Y0(1,63:93), num2cell(1:31), 'fontsize', 12, 'fontweight', 'bold', 'color', 'r')
% Yanim = Y0(1:4:end, :);
% for i = 2:size(Yanim, 1)
%     stickModify(handle, Yanim(i, :), connect);
%     pause(0.005)
% end
subplot(122), colormap('hot'), imagesc(Lambda), title('Ground truth network'), colorbar, daspect('manual')
figure(3), imagesc(Cy)
 % xyzankurAnimCompareMultiple(Ytest0, {}, 96, {'Ytest'}); % Animate mocap time-series.

%% Standard Glasso on mocap data, with varying lambda.
%{
warmLambda_hat = eye(d);
funObj = @(x)sparsePrecisionObj(x, d, nonZero, Cy);
for (i = 1:length(lambda))
    Lambda_hat{i} = eye(d);
    Lambda_hat{i}(nonZero) = L1GeneralProjection(funObj, warmLambda_hat(nonZero), lambda(i)*ones(d*d,1), options);
    GLASSOrocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat{i}<0);   % Performance evaluation.
    
    figure(3), imagesc(Lambda_hat{i}), colormap(hot), colorbar, title([ 'GLasso-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
end
TPs = GLASSOrocstats(:,1); FPs = GLASSOrocstats(:,2); FNs = GLASSOrocstats(:,3); TNs = GLASSOrocstats(:,4); 
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);

figure(2), clf, hold on, plot(Recalls, Precisions, '-xb'), text(Recalls, Precisions, num2cell(lambda)), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), legend('Glasso'), title('Recall-Precision')
figure(4), clf, hold on, plot(FPRs, Recalls, '-xb'), xlim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'GL auc: ' num2str(AUC) ], 4), title('ROC');
%}

%% Recovery of sparse-inverse and low-rank covariance via iterative application of GLASSO and RCA.
sigma2_n = 0.001 * trace(Cy);                                % Noise variance.
[S D] = eig(Cy);     [D perm] = sort(diag(D),'descend');    % Initialise W with a PPCA low-rank estimate.
W_hat_old = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));
WWt_hat_old = W_hat_old * W_hat_old';
Lambda_hat_old = eye(d); % pdinv(Cy); %
tic
parfor (i = 1:length(lambda),8)                             % Try different magnitudes of lambda.
    [WWt_hat_new, Lambda_hat_new, Lambda_hat_new_inv] = emrca(Y, WWt_hat_old, Lambda_hat_old, sigma2_n, lambda(i), nonZero, emrca_options);
    % Plot results.
    figure(5), clf, colormap('hot')
    subplot(131), imagesc(Lambda_hat_new), colorbar, title([ 'EM/RCA-recovered \Lambda with \lambda=', num2str(lambda(i)) ]);
    subplot(132), imagesc(Lambda_hat_new_inv), colorbar, title('\Sigma_{hat}'), colorbar
    subplot(133), imagesc(WWt_hat_new), colorbar, title('RCA-recovered WW'''), colorbar
    % Performance stats. Row format in pstats : [ TP FP FN TN ].
    EMRCArocstats(i,:) = emrcaRocStats(Lambda, Lambda_hat_new<0);
end
toc

%% Processing performance measures.
TPs = EMRCArocstats(:,1); FPs = EMRCArocstats(:,2); FNs = EMRCArocstats(:,3); TNs = EMRCArocstats(:,4); 
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
FPRs = FPs ./ (FPs + TNs);      AUC = trapz(flipud(FPRs), flipud(Recalls)) / max(FPRs);

%% Plot performance.
figure(3), plot(Recalls, Precisions, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision'), title('Recall-Precision')
hold on, text(Recalls, Precisions, num2cell(lambda)), legend('EM-RCA')
figure(4), hold on, plot(FPRs, Recalls, '-rs'), xlim([0 1]), ylim([0 1]), xlabel('FPR'), ylabel('TPR'), legend([ 'EM/RCA auc: ' num2str(AUC) ], 4), title('ROC')

