figure(1), clf
subplot(211)
plot(stats.alpha, stats.RMSerror, 'go-', 'Color', [0 .5 0],...
    'LineWidth', 3, 'MarkerFaceColor', [0 .5 0], 'MarkerSize', 5)
legend 'Iterative RCA'
ylim([3 3.5])
xlabel('\alpha', 'fontName', fontName, 'fontSize', fontSize)
ylabel('RMS error', 'fontName', fontName, 'fontSize', fontSize)
set(gca, 'fontname', fontName, 'fontsize', fontSize)

subplot(212)
plot(RMSerror_pCCA, 'ko-', 'LineWidth', 3, 'MarkerFaceColor',...
    [0 0 0], 'MarkerSize', 5)
legend('PCCA','location','southeast')
ylim([3 3.5])
xlabel('d', 'fontName', fontName, 'fontSize', fontSize)
ylabel('RMS error', 'fontName', fontName, 'fontSize', fontSize)
set(gca, 'fontname', fontName, 'fontsize', fontSize)

print -depsc ../tex/diagrams/RCAvsPCCA_noise3e-2.eps
system('epstopdf ../tex/diagrams/RCAvsPCCA_noise3e-2.eps')

figure(2), clf
% subplot(221),
hold on
plot(stats.alpha, stats.dz1, 'ro:', 'Color', [.75 0 0], 'LineWidth', 4, 'MarkerSize', 4)
plot(stats.alpha, stats.dz2, 'go-', 'Color', [0 .5 0], 'LineWidth', 4, 'MarkerSize', 4)
plot(stats.alpha, stats.dz3, 'bo--', 'Color', [0 0 .75], 'LineWidth', 4, 'MarkerSize', 4)
daspect([1 100 1])
xlabel('\alpha', 'fontName', fontName, 'fontSize', fontSize)
ylabel('Retained dimensions', 'fontName', fontName, 'fontSize', fontSize)
legend 'X1' 'Z' 'X2'
set(gca, 'fontname', fontName, 'fontsize', fontSize)

print -depsc ../tex/diagrams/RCAretdim_noise3e-2.eps
system('epstopdf ../tex/diagrams/RCAretdim_noise3e-2.eps')


% subplot(222), hold on,
% plot3(stats.alpha, stats.RMSerror, stats.dz1, 'r-'), grid on
% plot3(stats.alpha, stats.RMSerror, stats.dz2, 'g-')
% plot3(stats.alpha, stats.RMSerror, stats.dz3, 'b-')
% plot3(zeros(1,length(RMSerror_pCCA)), RMSerror_pCCA, 1:length(RMSerror_pCCA), 'k-')
% ylim([3 max(stats.RMSerror)])
% xlabel '\alpha', ylabel 'RMS error', zlabel 'Retained dimensions'
% legend 'z1' 'z2' 'z3' 'pCCA'
% pbaspect([1 1 1]), view([70 50]);

%{ 
save('demRCAsilhouette1DATA_NOISE03_2','stats')
cd ~/mlprojects/rca/tex/NIPS/diagrams/matlafigs/
saveas(gcf,'RCAvsPCCA_onSILHdata_CCAinit_Noise03_2', 'fig')
cd ~/mlprojects/rca/matlab/
printPlot('RCAvsPCCA_onSILHdata_CCAinit_noise03_2','../tex/NIPS/diagrams/')
%}


