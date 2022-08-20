function [ simulatedSpectra ] = simulate_spectra( targetSa, selectionParams, seedValue, nTrials )
% simulate response spectra from the target mean and covariance matrix

% Set initial seed for simulation
if seedValue ~= 0
    rng(seedValue);
else
    rng('default');
end

% Generate simulated response spectra with best matches to the target values
devTotalSim = zeros(nTrials,1);
for j=1:nTrials
    spectraSample{j} = exp(lhsnorm(targetSa.meanReq,targetSa.covReq,selectionParams.nGM));

    % evaluate simulation
    sampleMeanErr = mean(log(spectraSample{j})) - targetSa.meanReq; % how close is the mean of the spectra to the target
    sampleStdErr = std(log(spectraSample{j})) - sqrt(diag(targetSa.covReq))'; % how close is the standard dev. of the spectra to the target
    sampleSkewnessErr = skewness(log(spectraSample{j}),1); % how close is the skewness of the spectra to zero (i.e., the target)
    devTotalSim(j) = selectionParams.weights(1) * sum(sampleMeanErr.^2) + ...
                     selectionParams.weights(2) * sum(sampleStdErr.^2)+ ...
                     selectionParams.weights(3) * sum(sampleSkewnessErr.^2); % combine the three error metrics to compute a total error
end

[~, bestSample] = min(devTotalSim); % find the simulated spectra that best match the targets
simulatedSpectra = spectraSample{bestSample}; % return the best set of simulations

end
