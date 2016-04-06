function [ devTotal ] = compute_spectrum_error( selectionParams, targetSa, sampleSmall )
% compute the error between a selected set of spectra and the target
% spectrum



% Calculate the deviation from the target
if selectionParams.optType == 0
    if selectionParams.cond == 1 || (selectionParams.cond == 0 && selectionParams.isScaled == 0)
        % Compute deviations from target
        devMean = mean(sampleSmall) - targetSa.meanReq;
        devSig = std(sampleSmall) - targetSa.stdevs;
        devTotal = selectionParams.weights(1) * sum(devMean.^2) + selectionParams.weights(2) * sum(devSig.^2);
    end
    
    % Penalize bad spectra (set penalty to zero if this is not required)
    if selectionParams.penalty ~= 0
        for m=1:size(sampleSmall,2) % for each period
            devTotal = devTotal + sum(abs(exp(sampleSmall(:,m))>exp(targetSa.meanReq(m)+3*targetSa.stdevs(m)))) * selectionParams.penalty;
        end
    end
    
elseif selectionParams.optType == 1
    sortedlnSa = [min(sampleSmall); sort(sampleSmall)];
    norm_cdf = normcdf(sortedlnSa,repmat(targetSa.meanReq,selectionParams.nGM+1,1),repmat(targetSa.stdevs,selectionParams.nGM+1,1));
    emp_cdf = linspace(0,1,selectionParams.nGM+1);
    Dn = max(abs(repmat(emp_cdf',1,length(selectionParams.TgtPer)) - norm_cdf));
    devTotal = sum(Dn);
end




end

