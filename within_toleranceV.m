function [ withinTol, IMs ] = within_toleranceV( IMs, targetSa, selectionParams)
% check the mean and variance of the selected spectra versus their targets, 
% and see whether they are within tolerance

% Relative importance of components of GM
wH = 1-selectionParams.weightV;

% Extract data for H and V components of GM
nPer = length(selectionParams.TgtPer);
tgtSaMean = targetSa.meanReq(1,1:nPer);
tgtSaMeanV = targetSa.meanReq(1,(nPer+1):end);
tgtSaStd = targetSa.stdevs(1,1:nPer);
tgtSaStdV = targetSa.stdevs(1,(nPer+1):end);
selectedSa = IMs.sampleSmall;
selectedSaV = IMs.sampleSmallV;
stdevs = std(selectedSa);
stdevsV = std(selectedSaV);

% max error in median spectra
medianErrH = max(abs(exp(mean(selectedSa))- exp(tgtSaMean))./exp(tgtSaMean))*100;
medianErrV = max(abs(exp(mean(selectedSaV))- exp(tgtSaMeanV))./exp(tgtSaMeanV))*100;
medianErr = medianErrH*wH + medianErrV*(1-wH);
disp(['Max (across periods and components) error in median = ' num2str(medianErr,2) ' percent']);

% max error in standard deviation of log spectra
idNoTcond = ~ismember(selectionParams.TgtPer,selectionParams.Tcond);
stdErrH = max(abs(stdevs(idNoTcond) - tgtSaStd(idNoTcond))./tgtSaStd(idNoTcond))*100;
stdErrV = max(abs(stdevsV - tgtSaStdV)./tgtSaStdV)*100;
stdErr = stdErrH*wH + stdErrV*(1-wH);
disp(['Max (across periods and components) error in standard deviation = ' num2str(stdErr,2) ' percent']);

% Check whether errors are within the tolerance
withinTol = (medianErr < selectionParams.tol && stdErr < selectionParams.tol);
if withinTol
    disp(['The errors are within the target ' num2str(selectionParams.tol,2) ' percent tolerance']);
end

% Save errors for display
IMs.medianErr = medianErr;
IMs.stdErr = stdErr;

end

