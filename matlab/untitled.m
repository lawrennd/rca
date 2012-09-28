%%
for iSubject = 1:length(danceChannels)
    for iTrial = 1:length(danceChannels{iSubject})
        disp([iSubject, iTrial])
        skelPlayData(danceSkeletons{iSubject}{iTrial}, danceChannels{iSubject}{iTrial}, 1/120, Lambda_hat_emrca{5})
    end
end

cd ~/Desktop/CMUmocap/all_asfamc/05/
acclaimPlayFile('05.asf','05_08.amc')

skelPlayData(walkSkeletons{6}{1},walkChannels{6}{1}, 1/120, Lambda_hat_emrca{5}<0)      % Subject 9,    Trial 12
skelPlayData(runSkeletons{5}{1},runChannels{5}{1}, 1/120, Lambda_hat_emrca{5}<0)        % Subject 38,   Trial 3
skelPlayData(jumpSkeletons{3}{1},jumpChannels{3}{1}, 1/120, Lambda_hat_emrca{5}<0)      % Subject 49,   Trial 2
skelPlayData(danceSkeletons{1}{7}, danceChannels{1}{7}, 1/120, Lambda_hat_emrca{5}<0)   % Subject 5,    Trial 8,    rond de jambe in the air
skelPlayData(danceSkeletons{1}{4}, danceChannels{1}{4}, 1/120, Lambda_hat_emrca{5}<0)   % Subject 5,    Trial 5,    quasi-cou-de-pied
skelPlayData(danceSkeletons{3}{2}, danceChannels{3}{2}, 1/120, Lambda_hat_emrca{5}<0)   % Subject 55,   Trial 2,    lambada
skelPlayData(danceSkeletons{1}{6}, danceChannels{1}{6}, 1/120, Lambda_hat_emrca{5}<0)   % Subject 5,    Trial 7,    small-jetes


%%
nSamples = 10;
load DellaGattaData.mat
t_star = linspace(0,240,240)';
kern = kernCreate(t_star, 'rbf');
typDiff = 25; % Typical time difference between samples.
kern.inverseWidth = 1/(typDiff*typDiff); % Characteristic inverse-width.
% Covariance of the prior.
K_starStar = kernCompute(kern, t_star, t_star) + eye(length(t_star))*1e-6;
% Generate "training data" from a sine wave.
% t1 = rand(5, 1)*240;  t2 = rand(5, 1)*240;
f1 = sin(t1);    f2 = f1;
f2(end-1:end) = f2(end-1:end) + 3; 
f2([1 3]) = f([1 3]) - .4;

% Compute kernel for training data.
K_starf = kernCompute(kern, t_star, t1);
K_ff = kernCompute(kern, t1);
% Mean and covariance of posterior.
fbar = K_starf*pdinv(K_ff)*f1;
Sigma = K_starStar - K_starf*pdinv(K_ff)*K_starf';
% Sample from the posterior.
fsamp1 = real(gsamp(fbar, Sigma, nSamples));

% Compute kernel for training data.
K_starf = kernCompute(kern, t_star, t1);
K_ff = kernCompute(kern, t1);
% Mean and covariance of posterior.
fbar = K_starf*pdinv(K_ff)*f2;
Sigma = K_starStar - K_starf*pdinv(K_ff)*K_starf';
% Sample from the posterior.
fsamp2 = real(gsamp(fbar, Sigma, nSamples));

% Plot.
figure(2), clf, linHand = plot(t_star, fsamp1, 'k'); set(linHand, 'linewidth', 1)
% zeroAxes(gca, 0.025, 18, 'times')
hold on, plot(t1,f1, 'xb', 'markersize', 15, 'linewidth',3)
hold on, linHand = plot(t_star, fsamp2, 'k'); set(linHand, 'linewidth', 1)
% zeroAxes(gca, 0.025, 18, 'times')
plot(t1,f2, 'xr', 'markersize', 15, 'linewidth',3)
xlabel('time(m)'), ylabel('Gene Expression')

%% PCA on timeseries.
F = [fsamp1; fsamp2];
F = F - repmat(mean(F),size(F,1),1);
[V,D] = eig(F*F'./ size(F,1));
[D,ix] = sort(diag(D),'descend'); D=diag(D);
V = V(:,ix);
U = F'*V*sqrt(inv(D));
figure(5), plot(t_star,U'), legend(num2str(diag(D)))


%%
s = RandStream('mcg16807','Seed', 1e+6); RandStream.setDefaultStream(s)
K_starStar = kernCompute(kern, t_star, t_star) + eye(length(t_star))*1e-6;
fsamp = real(gsamp(zeros(size(t_star)), K_starStar, 2));
figure(2), clf, linHand = plot(t_star, fsamp); set(linHand, 'linewidth', 2)
zeroAxes(gca, 0.025, 18, 'times')

s = RandStream('mcg16807','Seed', 1e+6); RandStream.setDefaultStream(s)
X = randn(size(K_starStar),3)*.01;
fsamp = real(gsamp(zeros(size(t_star)), K_starStar + X*X', 2));
hold on, linHand = plot(t_star, fsamp); set(linHand, 'linewidth', 2)
zeroAxes(gca, 0.025, 18, 'times')
