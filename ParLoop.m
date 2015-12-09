function [ devTotal ] = ParLoop( devTotal, scaleFac, optInputs, sampleSmall, sampleBig, meanReq, stdevs, emp_cdf )
% Parallel loop to use within greedy optimization
optType = optInputs.optType;
if all(devTotal) && optType == 0 
    return;
end

PerTgt = optInputs.PerTgt;
cond = optInputs.cond;
isScaled = optInputs.isScaled;
weights = optInputs.weights;
penalty = optInputs.penalty;
maxScale = optInputs.maxScale;
recID = optInputs.recID;



parfor j = 1:optInputs.nBig
    sampleSmallTemp = [sampleSmall;sampleBig(j,:)+log(scaleFac(j))];
    
    % Calculate the appropriate measure of deviation and store in
    % devTotal (the D-statistic or the combination of mean and
    % sigma deviations)
    if optType == 0
        if cond == 1 || (cond == 0 && isScaled == 0)
            % Compute deviations from target
            devMean = mean(sampleSmallTemp) - meanReq;
            devSig = std(sampleSmallTemp) - stdevs;
            devTotal(j) = weights(1) * sum(devMean.^2) + weights(2) * sum(devSig.^2);
        end
        % Penalize bad spectra (set penalty to zero if this is not required)
        if penalty ~= 0
            for m=1:size(sampleSmall,1)
                devTotal(j) = devTotal(j) + sum(abs(exp(sampleSmallTemp(m,:))>exp(meanReq+3*stdevs'))) * penalty;
            end
        end
        
    elseif optType == 1
        [devTotal(j)] = KS_stat(PerTgt, emp_cdf, sampleSmallTemp, meanReq, stdevs);
    end
    
    % Scale factors for either type of optimization should not
    % exceed the maximum
    if (scaleFac(j) > maxScale)
        devTotal(j) = devTotal(j) + 1000000;
    end
    
    % Should cause improvement and record should not be repeated
    if (any(recID == j))
        devTotal(j) = 100000;
    end

end

end

