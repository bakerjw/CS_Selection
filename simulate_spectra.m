function [ simulatedSpectra ] = simulate_spectra( seedValue, nTrials, Tgts, optInputs )
% simulate response spectra from the target mean and covariance matrix

% Set initial seed for simulation
if seedValue ~= 0
    rng(seedValue); 
else
    rng('shuffle');
end


% find covariance values of zero and set them to a small number so that
% random number generation can be performed
Tgts.covReq( abs(Tgts.covReq) < 1e-17) = 1e-17;


% Generate simulated response spectra with best matches to the target values
devTotalSim = zeros(nTrials,1);
for j=1:nTrials
    spectraSample{j} = exp(lhsnorm(Tgts.meanReq,Tgts.covReq,optInputs.nGM));
    
    % evaluate simulation
    sampleMeanErr = mean(log(spectraSample{j})) - Tgts.meanReq; % how close is the mean of the spectra to the target
    sampleStdErr = std(log(spectraSample{j})) - sqrt(diag(Tgts.covReq))'; % how close is the standard dev. of the spectra to the target
    sampleSkewnessErr = skewness(log(spectraSample{j}),1); % how close is the skewness of the spectra to zero (i.e., the target)
    devTotalSim(j) = optInputs.weights(1) * sum(sampleMeanErr.^2) + ...
                     optInputs.weights(2) * sum(sampleStdErr.^2)+ ...
                     optInputs.weights(3) * sum(sampleSkewnessErr.^2); % combine the three error metrics to compute a total error
end

[~, bestSample] = min(devTotalSim); % find the simulated spectra that best match the targets 
simulatedSpectra = spectraSample{bestSample}; % return the best set of simulations

end

