% Diagnostic: Checking the RMS error and LML of the model so far.

C = blkdiag(W4W4,W1W1) + [W3; W2]*[W3' W2'] + Noise;     % Restriced RCA covariance model.
LML_RCA = -.5 * ((sy+sx)*log(2*pi) + logdet(C) + trace(sum(sum(pdinv(C)'.*S))));

YPred_RCA = (W2*W3'*pdinv(W4W4 + W3W3 + Nx)*XTest')'+ repmat(mean(Y,1), size(XTest,1), 1);
Ysqdiff = (YTest - YPred_RCA).^2;   RMSerror_RCA = sqrt(sum(Ysqdiff(:))/numel(Ysqdiff));

% subplot(221), imshow(C), title('C'), subplot(222), imshow(S), title('S')
% subplot(311), hold on, title('Trace approximation'), line([0 t],[trace(S) trace(S)]), tracesC(t)=trace(C); plot(1:t, tracesC, 'r-o'),  t=t+1;
% subplot(312), title('RMS error'), hold on,  subplot(313), title('LML'), hold on

disp(['####    Check trace balance: traceC = ' num2str(trace(C)) '    traceS = ' num2str(trace(S))])
disp(['####    LML_pCCA = ' num2str(LML_RCA)])
disp(['####    RMSerror_RCA = ' num2str(RMSerror_RCA)]);
disp(['####    dz1 = ' num2str(dz1) '    dz2 = ' num2str(dz2) '    dz3 = ' num2str(dz3)]);
