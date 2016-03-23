function [ withinTol ] = within_tolerance( selectedSa, targetSa, selectionParams)
% check the mean and variance of the selected spectra versus their targets, 
% and see whether they are within tolerance


% max error in mean log spectra
meanErr = max(abs(exp(mean(selectedSa))- exp(targetSa.meanReq))./exp(targetSa.meanReq))*100;
fprintf('Max (across periods) error in median = %3.1f percent \n', meanErr);

% max error in standard deviation of log spectra
stdevs = std(selectedSa);
stdErr = max(abs(stdevs(1:end ~= selectionParams.indTcond) - targetSa.stdevs(1:end ~= selectionParams.indTcond))./targetSa.stdevs(1:end ~= selectionParams.indTcond))*100;
fprintf('Max (across periods) error in standard deviation = %3.1f percent \n \n', stdErr);

% Check whether errors are within the tolerance
withinTol = (meanErr < selectionParams.tol && stdErr < selectionParams.tol);

if withinTol
    fprintf('The errors are within the required %3.1f percent tolerance \n', selectionParams.tol);
end

end

