% STABILITY Stability selection for Glasso and EM-RCA.
%
% FORMAT
% DESC
%
% SEEALSO : emrca
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2012
%
% RCA

clc, clear
s = RandStream('mcg16807','Seed', 1e+6); RandStream.setDefaultStream(s) % laplace 1e+6
nRestarts = 100;
stabilityThreshold = 0.5;
selection = randperm(2666);

%
results = cell(100,1);
for iRestart = 1:nRestarts
    disp(['Restart #', num2str(iRestart)]);
    
    demToyEMRca1                          % Insert demo script here for which you want to stabilise results.
    
    result = [];
%     result.Lambda_hat_glasso = Lambda_hat_glasso;
%     result.GLASSOrocstats = GLASSOrocstats;
    result.Lambda_hat_emrca  = Lambda_hat_emrca;
    result.EMRCArocstats = EMRCArocstats;
    results{iRestart} = result;
end
%}

%%
% stableLambda_idealglasso = cell(nReg_glasso, 1);
% IDEALGLASSOrocstats = zeros(nReg_glasso,4);
nReg_glasso = length(results{1}.Lambda_hat_glasso);
stableLambda_glasso = cell(nReg_glasso, 1);
GLASSOrocstats = zeros(nReg_glasso,4);
nReg_emrca = length(results{1}.Lambda_hat_emrca);
stableLambda_emrca = cell(nReg_emrca, 1);
EMRCArocstats = zeros(nReg_emrca,4);
%%
%
for iReg = 1:nReg_glasso
    stableLambda_glasso{iReg} = zeros(d);
    % stableLambda_idealglasso{iReg} = zeros(d);
    for iRestart = 1:nRestarts
        triuLambda_hat = triu(results{iRestart}.Lambda_hat_glasso{iReg,1}, 1);
        stableLambda_glasso{iReg} = stableLambda_glasso{iReg} + ...
            ... % (abs(triuLambda_hat) > 0);  % Edges in the estimated Lambda from Glasso.
            (triuLambda_hat < 0);
        
        % triuLambda_hat = triu(results{iRestart}.Lambda_hat_glasso{iReg,2}, 1);
        % stableLambda_idealglasso{iReg} = stableLambda_idealglasso{iReg} + ...
        %     (abs(triuLambda_hat) > 0);  % Edges in the estimated Lambda from Glasso, under non-confounded conditions.
    end
    stableLambda_glasso{iReg} = (stableLambda_glasso{iReg} >= round(nRestarts*stabilityThreshold));
    % stableLambda_idealglasso{iReg} = stableLambda_idealglasso{iReg} >= round(nRestarts*stabilityThreshold);
    GLASSOrocstats(iReg,:) = emrcaRocStats(Lambda, stableLambda_glasso{iReg});
    % IDEALGLASSOrocstats(iReg,:) = emrcaRocStats(Lambda, stableLambda_idealglasso{iReg});
    %     figure(3), imagesc(stableLambda_glasso{iReg}), colorbar
    %     figure(5), imagesc(stableLambda_idealglasso{iReg}), colorbar
end
TPs = GLASSOrocstats(:,1); FPs = GLASSOrocstats(:,2); FNs = GLASSOrocstats(:,3);
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
figure(2), hold on, plot(Recalls, Precisions, '--xb'), xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
% TPs = IDEALGLASSOrocstats(:,1); FPs = IDEALGLASSOrocstats(:,2); FNs = IDEALGLASSOrocstats(:,3);Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
% plot(Recalls, Precisions, '--om'), legend('GL','Ideal GL');
%}

%%
for iReg = 1:nReg_emrca
    stableLambda_emrca{iReg} = zeros(d);
    for iRestart = 1:nRestarts
        triuLambda_hat = triu(results{iRestart}.Lambda_hat_emrca{iReg}, 1);
        stableLambda_emrca{iReg} = stableLambda_emrca{iReg} + ...
            (abs(triuLambda_hat) > 0);          % Edges in the estimated Lambda from EM/RCA.
%             (triuLambda_hat < 0);             % Use for mocap demos.
    end
	stableLambda_emrca{iReg} = (stableLambda_emrca{iReg} >= round(nRestarts*stabilityThreshold));
    EMRCArocstats(iReg,:) = emrcaRocStats(Lambda, stableLambda_emrca{iReg});
end
TPs = EMRCArocstats(:,1); FPs = EMRCArocstats(:,2); FNs = EMRCArocstats(:,3);
Recalls = TPs ./ (TPs + FNs);   Precisions = TPs ./ (TPs + FPs);
figure(2), plot(Recalls, Precisions, '-rs'), hold on
% plot([1,.87,.27],[.275,.26,.57], 'g-', [1,.67,.2], [.275,.3,.5], 'b--')   % Literature performance.
xlim([0 1]), ylim([0 1]), xlabel('Recall'), ylabel('Precision')
% legend('EM-RCA','Kronecker-Glasso (reported)','Glasso (reported)')
legend('EM-RCA')

