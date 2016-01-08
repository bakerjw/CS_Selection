function [ sampleSmall, finalRecords, finalScaleFactors ] = GreedyOptPar( optInputs, Tgts, IMs )
% Parallelized greedy optimization

sampleSmall = IMs.sampleSmall;

display('Please wait...This algorithm takes a few minutes depending on the number of records to be selected');
if optInputs.cond == 0
    display('The algorithm is slower when scaling is used');
end
if optInputs.optType == 1
    display('The algorithm is slower when optimizing with the KS-test Dn statistic');
end

% if optimizing the ground motions by calculating the Dn value, first
% calculate the emperical CDF values (which will be the same at each
% period) and initialize a vector of Dn values
emp_cdf = 0;
if optInputs.optType == 1
    emp_cdf = linspace(0,1,optInputs.nGM+1);
end

numWorkers = 2;
parobj = parpool(numWorkers);

% Initialize scale factor vector
scaleFac = ones(optInputs.nBig,1);
for k=1:optInputs.nLoop % Number of passes
    
    for i=1:optInputs.nGM % Selects nGM ground motions
        display([num2str(round(((k-1)*optInputs.nGM + i-1)/(optInputs.nLoop*optInputs.nGM)*100)) '% done']);
        
        devTotal = zeros(optInputs.nBig,1);
        sampleSmall(i,:) = [];
        optInputs.recID(i,:) = [];
        
        if optInputs.isScaled == 1
            if optInputs.cond == 1
                scaleFac = exp(optInputs.lnSa1)./exp(IMs.sampleBig(:,optInputs.rec));
            elseif optInputs.cond == 0
                [scaleFac, devTotal] = bestScaleFactorPar(IMs.sampleBig, sampleSmall, Tgts.meanReq, Tgts.stdevs, optInputs.weights, optInputs.maxScale);
            end
        end
        
        % Try to add a new spectra to the subset list
        % new function 
        [devTotal] = ParLoop(devTotal, scaleFac, optInputs, sampleSmall, IMs.sampleBig, Tgts.meanReq,...
                             Tgts.stdevs, emp_cdf);                                
                                
        [minDevFinal, minID] = min(devTotal);
        % Add new element in the right slot
        if optInputs.isScaled == 1
            finalScaleFac(i) = scaleFac(minID);
        else
            finalScaleFac(i) = 1;
        end
        sampleSmall = [sampleSmall(1:i-1,:);IMs.sampleBig(minID,:)+log(scaleFac(minID));sampleSmall(i:end,:)];
        optInputs.recID = [optInputs.recID(1:i-1);minID;optInputs.recID(i:end)];
        
    end
    
    % Can the optimization be stopped after this loop based on the user
    % specified tolerance? Recalculate new standard deviations of new
    % sampleSmall and then recalculate new maximum percent errors of means
    % and standard deviations 
    if optInputs.optType == 0
        notT1 = find(optInputs.PerTgt ~= optInputs.PerTgt(optInputs.rec));
        stdevs = std(sampleSmall);
        meanErr = max(abs(exp(mean(sampleSmall))-Tgts.means)./Tgts.means)*100;
        stdErr = max(abs(stdevs(notT1) - Tgts.stdevs(notT1))./Tgts.stdevs(notT1))*100;
        fprintf('Max (across periods) error in median = %3.1f percent \n', meanErr); 
        fprintf('Max (across periods) error in standard deviation = %3.1f percent \n \n', stdErr);
        
        % If error is now within the tolerance, break out of the
        % optimization
        if meanErr < optInputs.tol && stdErr < optInputs.tol
            display('The percent errors between chosen and target spectra are now within the required tolerances.');
            break;
        end
    end
    
    fprintf('End of loop %1.0f of %1.0f \n', k, optInputs.nLoop) 
end

display('100% done');

% Output information
finalRecords = optInputs.recID;
finalScaleFactors = finalScaleFac';

delete(parobj);
end

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

function [ sumDn ] = KS_stat( periods, emp_cdf, sampleSmall, means, stdevs )
% calculate sum of all KS-test statistics 

sortedlnSa = [min(sampleSmall); sort(sampleSmall)];
norm_cdf = normcdf(sortedlnSa,repmat(means,size(sampleSmall,1)+1,1),repmat(stdevs,size(sampleSmall,1)+1,1));
Dn = max(abs(repmat(emp_cdf',1,length(periods)) - norm_cdf));
sumDn = sum(Dn);

end


