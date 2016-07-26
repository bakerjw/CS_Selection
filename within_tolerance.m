function [ withinTol ] = within_tolerance( selectedSa, targetSa, selectionParams)
% check the mean and variance of the selected spectra versus their targets, 
% and see whether they are within tolerance


% max error in median spectra
medianErr = max(abs(exp(mean(selectedSa))- exp(targetSa.meanReq))./exp(targetSa.meanReq))*100;
display(['Max (across periods) error in median = ' num2str(medianErr,2) ' percent']); 

% max error in standard deviation of log spectra
stdevs = std(selectedSa);
stdErr = max(abs(stdevs(1:end ~= selectionParams.indTcond) - targetSa.stdevs(1:end ~= selectionParams.indTcond))./targetSa.stdevs(1:end ~= selectionParams.indTcond))*100;
display(['Max (across periods) error in standard deviation = ' num2str(stdErr,2) ' percent']); 

% Check whether errors are within the tolerance
withinTol = (medianErr < selectionParams.tol && stdErr < selectionParams.tol);

if withinTol
    display(['The errors are within the target ' num2str(selectionParams.tol,2) ' percent tolerance']);
end

end

