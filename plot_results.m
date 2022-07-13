function [  ] = plot_results( selectionParams, targetSa, IMs, simulatedSpectra, SaKnown, knownPer )
%% Produce a number of plots of target and selected response spectra 

% Variables used here

% SaKnown    : As before, it contains the response spectra of all the
%              available ground motions (N*P matrix) - N ground motions,
%              P periods
% simulatedSpectra  : a matrix of simulated response spectra defined
%              only at PerTgt
% sampleBig  : Same as SaKnown, but is only defined at PerTgt, the
%              periods at which the target response spectrum properties
%              are computed
% sampleSmall: The response spectra of the selected ground motions,
%              defined at TgtPer
% meanReq    : Target mean for the (log) response spectrum
% covReq     : Target covariance for the (log) response spectrum

figureFontSize = 18; % font size for figure text

% Plot simulated response spectra
figure
loglog(selectionParams.TgtPer, exp(targetSa.meanReq), '-r', 'linewidth', 3)
hold on
loglog(selectionParams.TgtPer, exp(targetSa.meanReq + 1.96*sqrt(diag(targetSa.covReq))'), '--r', 'linewidth', 3)
loglog(selectionParams.TgtPer, simulatedSpectra', 'k');
loglog(selectionParams.TgtPer, exp(targetSa.meanReq - 1.96*sqrt(diag(targetSa.covReq))'), '--r', 'linewidth', 3)
axis([min(selectionParams.TgtPer) max(selectionParams.TgtPer) 1e-2 5])
xlabel('T (s)')
ylabel('S_a (g)')
legend('Median','2.5 and 97.5 percentile', 'Simulated spectra')
title('Response spectra of simulated ground motion spectra')
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)

% Plot selected response spectra
figure
loglog(selectionParams.TgtPer, exp(targetSa.meanReq), 'b', 'linewidth', 3)
hold on
loglog(selectionParams.TgtPer, exp(targetSa.meanReq + 1.96*sqrt(diag(targetSa.covReq))'), '--b', 'linewidth', 3)
loglog(knownPer,SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2)),'k');
loglog(selectionParams.TgtPer, exp(targetSa.meanReq - 1.96*sqrt(diag(targetSa.covReq))'), '--b', 'linewidth', 3)
axis([min(selectionParams.TgtPer) max(selectionParams.TgtPer) 1e-2 5])
xlabel('T (s)');
ylabel('S_a (g)');
legend('Median','2.5 and 97.5 percentile','Selected ground motions');
title ('Response spectra of selected ground motions');
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)

% Target, initial, and finally selected medians
figure
loglog(selectionParams.TgtPer, exp(targetSa.meanReq),'k','linewidth',1)
hold on
loglog(knownPer, exp(IMs.stageOneMeans),'r*--', 'linewidth',1)
loglog(knownPer,exp(mean(log(SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2))))),'b--','linewidth',1)
axis([min(selectionParams.TgtPer) max(selectionParams.TgtPer) 1e-2 5])
xlabel('T (s)');
ylabel('Median S_a (g)');
legend('Target', 'Stage 1 selection', 'Final selection');
title('Median Sa (i.e., exp(log mean Sa))')
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)

% Target, initial, and finally selected standard deviations
figure
semilogx(selectionParams.TgtPer,targetSa.stdevs,'k','linewidth',1)
hold on
semilogx(knownPer, IMs.stageOneStdevs,'r*--','linewidth',1)
semilogx(knownPer,std(log(SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2)))),'b--','linewidth',1)
axis([min(selectionParams.TgtPer) max(selectionParams.TgtPer) 0 1])
xlabel('T (s)');
ylabel('Standard deviation of lnS_a');
legend('Target', 'Stage 1 selection','Final selection');
title('Logarithmic Sa standard deviations')
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)


% Compute sample correlations from selected spectra
selectedSa = log(SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2)));
selectedCorr = corrcoef(selectedSa);

% Calculate target correlations using the target covariance matrix
knownStdevs = sqrt(diag(targetSa.covAllT));
corrReq = targetSa.covAllT./(knownStdevs*knownStdevs');

% Contours of correlations from selected spectra
figure
contour(knownPer, knownPer, selectedCorr);
set(gca,'yscale','log','xscale','log');
axis square;
title('Sample correlation coefficients contour');
xlabel('T_1 (s)')
ylabel('T_2 (s)')
colorbar('YLim',[0 1]);
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)

% Contours of target correlations
figure
contour(knownPer, knownPer, corrReq);
set(gca,'yscale','log','xscale','log');
axis square;
xlabel('T_1 (s)');
ylabel('T_2 (s)');
title('Target correlation coefficients contour');
colorbar('YLim',[0 1]);
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)

% Difference between target and sample correlations
diffCorr = abs(selectedCorr-corrReq);
figure
contour(knownPer, knownPer,diffCorr);
set(gca, 'yscale', 'log','xscale','log');
axis square;
title('abs(sample correlation - target correlation)');
xlabel('T_1 (s)');
ylabel('T_2 (s)');
colorbar('YLim',[0 1]);
caxis([0 1]);
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)


