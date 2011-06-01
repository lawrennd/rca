clf
subplot(221), hold on
plot(stats.alpha, stats.dz1, 'r.-')
plot(stats.alpha, stats.dz2, 'g.-')
plot(stats.alpha, stats.dz3, 'b.-')
xlabel '\alpha', ylabel 'Retained dimensions', legend 'z1' 'z2' 'z3'
ylim([0 max([stats.dz1 stats.dz2 stats.dz3])]);

subplot(223), plot(stats.alpha, stats.RMSerror, 'g.-')
xlabel '\alpha', ylabel 'RMS error', legend 'RCA'
% ylim([min([stats.RMSerror RMSerror_pCCA]) max([stats.RMSerror RMSerror_pCCA])])
ylim([min([stats.RMSerror RMSerror_pCCA]) 3.5])

subplot(224), plot(RMSerror_pCCA, 'k.-')
xlabel 'dim', ylabel 'RMS error', legend 'pCCA'
% ylim([min([stats.RMSerror RMSerror_pCCA]) max([stats.RMSerror RMSerror_pCCA])])
ylim([min([stats.RMSerror RMSerror_pCCA]) 3.5])

subplot(222), hold on,
plot3(stats.alpha, stats.RMSerror, stats.dz1, 'r-'), grid on
plot3(stats.alpha, stats.RMSerror, stats.dz2, 'g-')
plot3(stats.alpha, stats.RMSerror, stats.dz3, 'b-')
plot3(zeros(1,length(RMSerror_pCCA)), RMSerror_pCCA, 1:length(RMSerror_pCCA), 'k-')
ylim([3 max(stats.RMSerror)])
xlabel '\alpha', ylabel 'RMS error', zlabel 'Retained dimensions'
legend 'z1' 'z2' 'z3' 'pCCA'
pbaspect([1 1 1]), view([70 50]);

%{ 
save('demRCAsilhouette1DATA_NOISE03_2','stats')
cd ~/mlprojects/rca/tex/NIPS/diagrams/matlafigs/
saveas(gcf,'RCAvsPCCA_onSILHdata_CCAinit_Noise03_2', 'fig')
cd ~/mlprojects/rca/matlab/
printPlot('RCAvsPCCA_onSILHdata_CCAinit_noise03_2','../tex/NIPS/diagrams/')
%}


