function [ devTotal, withinTol ] = compute_spectrum_error( selectionParams, targetSa, sampleSmall )
% compute the error between a selected set of spectra and the target
% spectrum



% Calculate the deviation from the target
if selectionParams.optType == 0
    if selectionParams.cond == 1 || (selectionParams.cond == 0 && selectionParams.isScaled == 0)
        % Compute deviations from target
        sampleMean = sum(sampleSmall)/size(sampleSmall,1);
        sampleStd = sqrt(sum((sampleSmall - ones(size(sampleSmall,1),1)*sampleMean).^2) / size(sampleSmall,1));
        devMean = abs(exp(sampleMean) - exp(targetSa.meanReq)) ./ exp(targetSa.meanReq);
        devStd =  abs(sampleStd  - targetSa.stdevs)  ./ targetSa.stdevs;
        devTotal = (1-selectionParams.stdWeight) * max(devMean.^2) + selectionParams.stdWeight * max(devStd.^2);
    end
    
    % Penalize bad spectra (if penalty has been set to non-zero)
    if selectionParams.penalty ~= 0
        for m=1:size(sampleSmall,2) % for each period
            devTotal = devTotal + sum(abs(exp(sampleSmall(:,m))>exp(targetSa.meanReq(m)+3*targetSa.stdevs(m)))) * selectionParams.penalty;
        end
    end
    
    % Check whether errors are within the tolerance
    withinTol = (devTotal < selectionParams.tol);

elseif selectionParams.optType == 1
    sortedlnSa = [min(sampleSmall); sort(sampleSmall)];
    norm_cdf = normcdf(sortedlnSa,repmat(targetSa.meanReq,selectionParams.nGM+1,1),repmat(targetSa.stdevs,selectionParams.nGM+1,1));
    emp_cdf = linspace(0,1,selectionParams.nGM+1);
    Dn = max(abs(repmat(emp_cdf',1,length(selectionParams.TgtPer)) - norm_cdf));
    devTotal = mean(Dn);
    
    % Check whether errors are within the tolerance
    withinTol = (devTotal < selectionParams.dStatTol);
end




end

