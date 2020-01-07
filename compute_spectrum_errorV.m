function [ devTotal ] = compute_spectrum_errorV( selectionParams, targetSa, testSpectra )
% compute the error between a selected set of spectra and the target
% spectrum

% Extract data for V component of GM
nPer = length(selectionParams.TgtPer);
sampleSmall = testSpectra(:,1:nPer);
sampleSmallV = testSpectra(:,(nPer+1):end);
tgtSaMean = targetSa.meanReq(1,1:nPer);
tgtSaMeanV = targetSa.meanReq(1,(nPer+1):end);
tgtSaStd = targetSa.stdevs(1,1:nPer);
tgtSaStdV = targetSa.stdevs(1,(nPer+1):end);

% Relative importance of components of GM
wH = 1-selectionParams.weightV;

% Calculate the deviation from the target
if selectionParams.optType == 0 % computed weighted sum of squared errors of mean and standard deviation
    % Compute deviations from target (H component)
    sampleMean = sum(sampleSmall)/size(sampleSmall,1);
    sampleVar = sum((sampleSmall - ones(size(sampleSmall,1),1)*sampleMean).^2) / size(sampleSmall,1);
    devMean = sampleMean - tgtSaMean;
    devSig = sqrt(sampleVar) - tgtSaStd;
    
    % Compute deviations from target (V component)
    sampleMeanV = sum(sampleSmallV)/size(sampleSmallV,1);
    sampleVarV = sum((sampleSmallV - ones(size(sampleSmallV,1),1)*sampleMeanV).^2) / size(sampleSmallV,1);
    devMeanV = sampleMeanV - tgtSaMeanV;
    devSigV = sqrt(sampleVarV) - tgtSaStdV;
    
    % Combine using wH
    devTotal = selectionParams.weights(1) * (sum(devMean.^2)*wH + sum(devMeanV.^2)*(1-wH)) + selectionParams.weights(2) * (sum(devSig.^2)*wH + sum(devSigV.^2)*(1-wH));
    
    % Penalize bad spectra (set penalty to zero if this is not required)
    if selectionParams.penalty ~= 0
        for m=1:size(sampleSmall,2) % for each period
            devTotal = devTotal + sum(abs(exp(sampleSmall(:,m))>exp(tgtSaMean(m)+3*tgtSaStd(m)))) * selectionParams.penalty;
        end
        for m=1:size(sampleSmallV,2) % repeat penalty for V component
            devTotal = devTotal + sum(abs(exp(sampleSmallV(:,m))>exp(tgtSaMeanV(m)+3*tgtSaStdV(m)))) * selectionParams.penalty;
        end
    end
elseif selectionParams.optType == 1 % compute Kolmogorov-Smirnov error metric
    emp_cdf = linspace(0,1,selectionParams.nGM+1);
    % KS metric for H component only
    sortedlnSa = [min(sampleSmall); sort(sampleSmall)];
    norm_cdf = normcdf(sortedlnSa,repmat(tgtSaMean,selectionParams.nGM+1,1),repmat(tgtSaStd,selectionParams.nGM+1,1));
    Dn = max(abs(repmat(emp_cdf',1,length(selectionParams.TgtPer)) - norm_cdf));
    % KS metric for V component only
    sortedlnSaV = [min(sampleSmallV); sort(sampleSmallV)];
    norm_cdfV = normcdf(sortedlnSaV,repmat(tgtSaMeanV,selectionParams.nGM+1,1),repmat(tgtSaStdV,selectionParams.nGM+1,1));
    DnV = max(abs(repmat(emp_cdf',1,length(selectionParams.TgtPerV)) - norm_cdfV));
    % Combine for H and V components    
    devTotal = sum(Dn)*wH + sum(DnV)*(1-wH);
end




end

