function [  ] = plot_resultsV( selectionParams, targetSa, IMs, simulatedSpectra, SaKnown, knownPer )
%% Produce a number of plots of target and selected response spectra 

% Variables used here

% SaKnown           : As before, it contains the response spectra of all the
%                   available ground motions (N*P matrix) - N ground motions,
%                   P periods
% SaKnownV          : Same as SaKnown except for V component
% simulatedSpectra  : a matrix of simulated response spectra defined
%                   at selectionParams.TgtPer and selectionParams.TgtPerV
% sampleBig         : Same as SaKnown, but is only defined at PerTgt, the
%                   periods at which the target response spectrum properties
%                   are computed
% sampleBigV        : Same as sampleBig except for V component
% sampleSmall       : The response spectra of the selected ground motions,
%                   defined at TgtPer
% sampleSmallV      : Same as sampleSmall except for V component
% meanReq    : Target mean for the (log) response spectrum
% covReq     : Target covariance for the (log) response spectrum

% Extract data for H and V components of GM
nPer = length(selectionParams.TgtPer);
tgtSpecMed = exp(targetSa.meanReq(1,1:nPer));
tgtSpecMedV = exp(targetSa.meanReq(1,(nPer+1):end));
tgtSpecStd = targetSa.stdevs(1,1:nPer);
tgtSpecStdV = targetSa.stdevs(1,(nPer+1):end);
simSpec = simulatedSpectra(:,1:nPer);
simSpecV = simulatedSpectra(:,(nPer+1):end);
SaKnownV = selectionParams.SaKnownV;


% Plot simulated response spectra
nSTD = 1.96;
figure
subplot(1,2,1);
loglog(selectionParams.TgtPer, tgtSpecMed, '-r', 'linewidth', 3)
hold on
loglog(selectionParams.TgtPer, tgtSpecMed .* exp(+nSTD*tgtSpecStd), '--r', 'linewidth', 3)
for i=1:selectionParams.nGM
    loglog(selectionParams.TgtPer, simSpec(i,:), 'k');
end
loglog(selectionParams.TgtPer, tgtSpecMed .* exp(-nSTD*tgtSpecStd), '--r', 'linewidth', 3)
axis([min(selectionParams.TgtPer) max(selectionParams.TgtPer) 1e-2 5])
axis square
xlabel('T (s)')
ylabel('S_a (g)')
hLeg = legend('Median','2.5 and 97.5 percentile', 'Simulated spectra');
set(hLeg,'location','southwest','box','off');
title('Horizontal')

subplot(1,2,2);
loglog(selectionParams.TgtPerV, tgtSpecMedV, '-r', 'linewidth', 3)
hold on
loglog(selectionParams.TgtPerV, tgtSpecMedV .* exp(+nSTD*tgtSpecStdV), '--r', 'linewidth', 3)
for i=1:selectionParams.nGM
    loglog(selectionParams.TgtPerV, simSpecV(i,:), 'k');
end
loglog(selectionParams.TgtPerV, tgtSpecMedV .* exp(-nSTD*tgtSpecStdV), '--r', 'linewidth', 3)
axis([min(selectionParams.TgtPerV) max(selectionParams.TgtPerV) 1e-2 5])
axis square
xlabel('T (s)')
title('Vertical')


% Plot selected response spectra
figure
subplot(1,2,1);
loglog(selectionParams.TgtPer, tgtSpecMed, 'b', 'linewidth', 3)
hold on
loglog(selectionParams.TgtPer, tgtSpecMed .* exp(+nSTD*tgtSpecStd), '--b', 'linewidth', 3)
loglog(knownPer,SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2)),'k');
loglog(selectionParams.TgtPer, tgtSpecMed .* exp(-nSTD*tgtSpecStd), '--b', 'linewidth', 3)
axis([min(selectionParams.TgtPer) max(selectionParams.TgtPer) 1e-2 5])
axis square
xlabel('T (s)');
ylabel('S_a (g)');
hLeg = legend('Median','2.5 and 97.5 percentile','Selected ground motions');
set(hLeg,'location','southwest','box','off');
title('Horizontal')

subplot(1,2,2);
loglog(selectionParams.TgtPerV, tgtSpecMedV, 'b', 'linewidth', 3)
hold on
loglog(selectionParams.TgtPerV, tgtSpecMedV .* exp(+nSTD*tgtSpecStdV), '--b', 'linewidth', 3)
loglog(knownPer,SaKnownV(IMs.recID,:).*repmat(IMs.scaleFacV,1,size(SaKnownV,2)),'k');
loglog(selectionParams.TgtPerV, tgtSpecMedV .* exp(-nSTD*tgtSpecStdV), '--b', 'linewidth', 3)
axis([min(selectionParams.TgtPerV) max(selectionParams.TgtPerV) 1e-2 5])
axis square
xlabel('T (s)');
title('Vertical')


% Target, initial, and finally selected medians
figure
subplot(1,2,1);
loglog(selectionParams.TgtPer, tgtSpecMed,'k','linewidth',1)
hold on
loglog(knownPer, exp(IMs.stageOneMeans),'r*--', 'linewidth',1)
loglog(knownPer,exp(mean(log(SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2))))),'b--','linewidth',1)
axis([min(selectionParams.TgtPer) max(selectionParams.TgtPer) 1e-2 5])
axis square
xlabel('T (s)');
ylabel('Median S_a (g)');
hLeg = legend('Target', 'Stage 1 selection', 'Final selection');
set(hLeg,'location','southwest','box','off');
title('Horizontal')

subplot(1,2,2);
loglog(selectionParams.TgtPerV, tgtSpecMedV,'k','linewidth',1)
hold on
loglog(knownPer, exp(IMs.stageOneMeansV),'r*--', 'linewidth',1)
loglog(knownPer,exp(mean(log(SaKnownV(IMs.recID,:).*repmat(IMs.scaleFacV,1,size(SaKnownV,2))))),'b--','linewidth',1)
axis([min(selectionParams.TgtPerV) max(selectionParams.TgtPerV) 1e-2 5])
axis square
xlabel('T (s)');
title('Vertical')


% Target, initial, and finally selected standard deviations
figure
subplot(1,2,1);
semilogx(selectionParams.TgtPer,tgtSpecStd,'k','linewidth',1)
hold on
semilogx(knownPer, IMs.stageOneStdevs,'r*--','linewidth',1)
semilogx(knownPer,std(log(SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2)))),'b--','linewidth',1)
axis([min(selectionParams.TgtPer) max(selectionParams.TgtPer) 0 1])
axis square
xlabel('T (s)');
ylabel('Standard deviation of lnS_a');
hLeg = legend('Target', 'Stage 1 selection','Final selection');
set(hLeg,'location','northeast','box','off');
title('Horizontal')

subplot(1,2,2);
semilogx(selectionParams.TgtPerV,tgtSpecStdV,'k','linewidth',1)
hold on
semilogx(knownPer, IMs.stageOneStdevsV,'r*--','linewidth',1)
semilogx(knownPer,std(log(SaKnownV(IMs.recID,:).*repmat(IMs.scaleFacV,1,size(SaKnownV,2)))),'b--','linewidth',1)
axis([min(selectionParams.TgtPerV) max(selectionParams.TgtPerV) 0 1])
axis square
xlabel('T (s)');
title('Vertical')