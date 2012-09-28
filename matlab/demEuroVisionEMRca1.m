% DEMEUROVISIONRCA1 EM-RCA demo on the recovery of the Eurovision collusion
% network.
%
% FORMAT
% DESC
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2012
%
% RCA

clc
addpath(genpath('~/mlprojects/matlab/general/'))
importTool({'rca','ndlutil'})
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
addpath(genpath('~/mlprojects/rca/matlab/EUROVISION_DATA'))
addpath(genpath('~/mlprojects/rca/matlab/graphViz4Matlab'))


% General settings.
colormap hot
options = struct('verbose',1,'order',-1);
allCountryNames = {'Albania', 'Andorra', 'Armenia', 'Austria', ...
    'Azerbaijan', 'Belarus', 'Belgium', 'Bosnia-Herzegovina', ...
    'Bulgaria', 'Croatia', 'Cyprus', 'Czech Republic', ...
    'Denmark', 'Estonia', 'Finland', 'France', ...
    'Georgia', 'Germany', 'Greece', 'Hungary', ...
    'Iceland', 'Ireland', 'Israel', 'Italy', ...
    'Latvia', 'Lithuania', 'Macedonia', 'Malta', ...
    'Moldova', 'Monaco', 'Montenegro', 'Netherlands', ...
    'Norway', 'Poland', 'Portugal', 'Romania', ...
    'Russia', 'San Marino', 'Serbia', 'Serbia and Montenegro', ...
    'Slovakia', 'Slovenia', 'Spain', 'Sweden', ...
    'Switzerland', 'Turkey', 'Ukraine', 'United Kingdom'};
allCountryCoordinates = ...
    fliplr([41,20;  42.55,1.61;  40.07,45.04;  47.51,14.55;
    40.48,48.01;  53.749,27.905;  50.695,4.68;  44.55,17.71;
    42.8,25.291;  45.198,15.359;  35.126,33.43;  49.85,15.315;
    55.974,10.1;  58.768,25.708;  64.09,26.543;  46.228,2.214;
    42.033,43.767;  51.166,10.45;  39.074,21.824;  47.16,19.5;
    64.96,-19;  53.413,-8.244;  31.046,34.852;  41.87,12.57;
    56.88,24.6;  55.17,23.88;  41.61,21.745;  35.94,14.375;
    47.41,28.37;  43.0384,7.4246;  42.709,19.374;  52.133,5.291;
    60.47,8.469;  51.92,19.145;  39.4,-8.224;  45.94,24.97;
    55.743,37.615;  43.94236,12.4577;  44.017,21;  43.1651,19.7095;
    48.67,19.7;  46.151,14.995;  40.464,-3.75;  60.13,18.644;
    46.8182,8.228;  38.96,35.24;  48.38,31.17;  55.378,-3.436]);
allCountryCoordinates = allCountryCoordinates - repmat( min(allCountryCoordinates), size(allCountryCoordinates,1),1 );


%% Import data.
[temp, Yfinals] = importEVfile('./EUROVISION_DATA/finals1998-2012.csv', allCountryNames, 'active');
[temp, Ysemifinals] = importEVfile('./EUROVISION_DATA/semifinals2004-2012.csv', allCountryNames, 'active');
Y = [Yfinals; Ysemifinals];

% Subtract mean and standardize.
Y = Y - repmat(mean(Y), size(Y,1), 1);
Y = Y ./ repmat(std(Y), size(Y,1), 1);
Cy = Y' * Y / size(Y,1);
[n,d] = size(Y);

%% Recovery of sparse-inverse and low-rank covariance via iterative application of GLASSO and RCA.
lambda = 5.^linspace(-8,-5,10);  % lambda = 5.^linspace(-8,3,30);
emrca_options = struct('limit',1e-4, 'showProgress',1 , 'verbose',0, 'errorCheck',0, 'maxNumIter',1000);
[S D] = eig(Cy);     [D perm] = sort(diag(D),'descend');    % Initialise W with a PCA low-rank estimate.
sigma2_n = D(2); % 3*trace(Cy)/d;
W_hat_old = S(:,perm(D>sigma2_n)) * sqrt(diag(D(D>sigma2_n)-sigma2_n));
% W_hat_old = S(:,perm(1)) * sqrt(D(1)-sigma2_n);
WWt_hat_old = W_hat_old * W_hat_old';
Lambda_hat_emrca = cell(length(lambda),1);
W_hat_new = cell(length(lambda),1);
Lambda_hat_old = eye(d);
nonZero = find(ones(d));                    % Induce any prior knowledge of zeros. If not, this is all ones.
tic
parfor (i = 1:length(lambda),8)            % Try different magnitudes of lambda.
    [W_hat_new{i}, Lambda_hat_emrca{i}] = emrca(Y, WWt_hat_old, Lambda_hat_old, sigma2_n, lambda(i), nonZero, emrca_options);
end
toc

for i = 1:length(lambda)
    % Plot results.
    figure(1), subplot(5,6,2*i-1), imagesc(Lambda_hat_emrca{i}~=0), title([ '\lambda=', num2str(lambda(i)) ]);
    subplot(5,6,2*i), imagesc(W_hat_new{i}, [-.5 .5]), title('W')
end

figure(3), clf, colormap hot, k=7; Lambda = zeros(d); Lambda(Lambda_hat_emrca{k}<0) = Lambda_hat_emrca{k}(Lambda_hat_emrca{k}<0);
% subplot(1,10,1:7)
visualiseNetwork(Lambda, 'edgeColor',Lambda, 'nodeLabels', allCountryNames, 'nodePositions',allCountryCoordinates, 'fontSize',6, 'nodeSize',0);
printPlot('EuroNet')
% [Ws, Wsi] = sort(W_hat_new{k},'descend');
% subplot(1,10,9:10), imagesc(Ws), set(gca, 'YTick',1:d, 'YTickLabel', allCountryNames(Wsi)), axis image,
% colorbar

% title([ 'Recovered \Lambda with \lambda=', num2str(lambda(9)) ]),
% set(gca, 'YTick', 1:48, 'YTickLabel', allCountryNames), xticklabel_rotate(1:48, 90, allCountryNames);

% g = graphViz4Matlab('-adjMat',Lambda_hat_emrca{5}, '-nodeLabels',allCountryNames, '-undirected',true);
% g.setNodePositions(allCountryCoordinates);
