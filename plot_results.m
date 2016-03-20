function [  ] = plot_results( optInputs, Tgts, IMs, simulatedSpectra, SaKnown, knownPer, knownCovReq )
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
loglog(optInputs.TgtPer, exp(Tgts.meanReq), '-r', 'linewidth', 3)
hold on
loglog(optInputs.TgtPer, exp(Tgts.meanReq + 1.96*sqrt(diag(Tgts.covReq))'), '--r', 'linewidth', 3)
loglog(optInputs.TgtPer, simulatedSpectra, 'k');
loglog(optInputs.TgtPer, exp(Tgts.meanReq - 1.96*sqrt(diag(Tgts.covReq))'), '--r', 'linewidth', 3)
axis([min(optInputs.TgtPer) max(optInputs.TgtPer) 1e-2 5])
xlabel('T (s)')
ylabel('S_a (g)')
legend('Median','2.5 and 97.5 percentile', 'Simulated spectra')
title('Response spectra of simulated ground motion spectra')
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)

% Plot selected response spectra
figure
loglog(optInputs.TgtPer, exp(Tgts.meanReq), 'b', 'linewidth', 3)
hold on
loglog(optInputs.TgtPer, exp(Tgts.meanReq + 1.96*sqrt(diag(Tgts.covReq))'), '--b', 'linewidth', 3)
loglog(knownPer,SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2)),'k');
loglog(optInputs.TgtPer, exp(Tgts.meanReq - 1.96*sqrt(diag(Tgts.covReq))'), '--b', 'linewidth', 3)
axis([min(optInputs.TgtPer) max(optInputs.TgtPer) 1e-2 5])
xlabel('T (s)');
ylabel('S_a (g)');
legend('Median','2.5 and 97.5 percentile','Selected ground motions');
title ('Response spectra of selected ground motions');
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)

% Target, initial, and finally selected medians
figure
loglog(optInputs.TgtPer, exp(Tgts.meanReq),'k','linewidth',1)
hold on
loglog(knownPer, exp(IMs.stageOneMeans),'r*--', 'linewidth',1)
loglog(knownPer,exp(mean(log(SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2))))),'b--','linewidth',1)
axis([min(optInputs.TgtPer) max(optInputs.TgtPer) 1e-2 5])
xlabel('T (s)');
ylabel('Median S_a (g)');
legend('Target', 'Stage 1 selection', 'Stage 2 selection');
title('Target and sample exponential logarithmic means (i.e., medians)')
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)

% Target, initial, and finally selected standard deviations
figure
semilogx(optInputs.TgtPer,Tgts.stdevs,'k','linewidth',1)
hold on
semilogx(knownPer, IMs.stageOneStdevs,'r*--','linewidth',1)
semilogx(knownPer,std(log(SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2)))),'b--','linewidth',1)
axis([min(optInputs.TgtPer) max(optInputs.TgtPer) 0 1])
xlabel('T (s)');
ylabel('Standard deviation of lnS_a');
legend('Target', 'Stage 1 selection','Final selection');
title('Target and sample logarithmic standard deviations')
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)


% Compute sample correlations from selected spectra
sampleUse = log(SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2)));
sampleUse = [sampleUse(:,knownPer<optInputs.T1) interp1(knownPer,sampleUse',optInputs.T1)' sampleUse(:,knownPer>optInputs.T1 & knownPer<=10)];

sampleCorr = corrcoef(sampleUse);

% sampleCorr = zeros(length(knownPer(knownPer<=10)));
% for i=1:length(knownPer(knownPer<=10))
%     for j=1:length(knownPer(knownPer<=10))
%         corrMatrix = corrcoef((sampleUse(:,i)),(sampleUse(:,j)));
%         sampleCorr(i,j) = corrMatrix(1,2);
%     end
% end



% Calculate target correlations from covariances
knownStdevs = sqrt(diag(knownCovReq));
corrReq = knownCovReq./(knownStdevs*knownStdevs');

% Contours of correlations from selected spectra
figure
contour(knownPer(knownPer<=10), knownPer(knownPer<=10), sampleCorr);
set(gca,'yscale','log','xscale','log');
axis square;
xlabel('T_1');
ylabel('T_2');
title('Sample correlation coefficients contour');
xlabel('T_1 (s)')
ylabel('T_2 (s)')
colorbar('YLim',[0 1]);
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)

% Contours of target correlations
figure
contour(knownPer(knownPer<=10), knownPer(knownPer<=10), corrReq(knownPer<=10, knownPer<=10));
set(gca,'yscale','log','xscale','log');
axis square;
xlabel('T_1 (s)');
ylabel('T_2 (s)');
title('Target correlation coefficients contour');
colorbar('YLim',[0 1]);
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)

% Difference between target and sample correlations
diffCorr = abs(sampleCorr-corrReq(knownPer<=10,knownPer<=10));
figure
contour(knownPer(knownPer<=10),knownPer(knownPer<=10),diffCorr);
set(gca, 'yscale', 'log','xscale','log');
axis square;
title('Difference in the correlation (sample-target)');
xlabel('T_1 (s)');
ylabel('T_2 (s)');
colorbar('YLim',[0 1]);
caxis([0 1]);
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)


