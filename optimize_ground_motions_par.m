function [ IMs ] = optimize_ground_motions_par( selectionParams, targetSa, IMs )
% Parallelized greedy optimization, for variable definitions, see
% optimize_ground_motions(selectionParams, targetSa, IMs)

% Note that launching the parallel pool takes some time on its own, so it
% is typically only worth using this parallel code if there are large 
% numbers of ground motions to select (i.e. over 100) 

% Note also that this function is somewhat experimental. It has not been
% optimized as completely as optimize_ground_motions.m, but it should
% produce equivalent results.

sampleSmall = IMs.sampleSmall;

if selectionParams.cond == 0 && selectionParams.isScaled
    display('The algorithm is slower when scaling is used');
end
if selectionParams.optType == 1
    display('The algorithm is slower when optimizing with the KS-test Dn statistic');
end

% if optimizing the ground motions by calculating the Dn value, first
% calculate the emperical CDF values (which will be the same at each
% period) and initialize a vector of Dn values
emp_cdf = 0;
if selectionParams.optType == 1
    emp_cdf = linspace(0,1,selectionParams.nGM+1);
end

numWorkers = 2; % specify the number of workers
parobj = parpool(numWorkers);

% Initialize scale factor vectors if possible
if selectionParams.isScaled == 0 % no scaling so set scale factors = 1
    scaleFac = ones(selectionParams.nBig,1);
    IMs.scaleFac = ones(selectionParams.nGM,1);
elseif selectionParams.isScaled && selectionParams.cond % Sa(Tcond) scaling
    scaleFac = exp(selectionParams.lnSa1)./exp(IMs.sampleBig(:,selectionParams.indTcond));
end
    
hw = waitbar(0,'Optimizing ground motion selection');

for k=1:selectionParams.nLoop % Number of passes
    for i=1:selectionParams.nGM % Selects nGM ground motions
        
        devTotal = zeros(selectionParams.nBig,1);
        sampleSmall(i,:) = [];
        IMs.recID(i,:) = [];
        
        % if scaling with unconditional selection, compute scale factors
        if selectionParams.isScaled && selectionParams.cond == 0
            scaleFac = compute_scale_factor(IMs.sampleBig, sampleSmall, targetSa.meanReq, targetSa.stdevs, selectionParams.weights, selectionParams.maxScale);
        end
        
        % Try to add a new spectra to the subset list
        % new function 
        [devTotal] = ParLoop(devTotal, scaleFac, selectionParams, sampleSmall, IMs, targetSa.meanReq,...
                             targetSa.stdevs, emp_cdf);                                
                                
        [~, minID] = min(devTotal);
        % Add new element in the right slot
        IMs.scaleFac(i) = scaleFac(minID);
        IMs.recID = [IMs.recID(1:i-1);minID;IMs.recID(i:end)];
        sampleSmall = [sampleSmall(1:i-1,:);IMs.sampleBig(minID,:)+log(scaleFac(minID));sampleSmall(i:end,:)];
        
        waitbar(((k-1)*selectionParams.nGM + i)/(selectionParams.nLoop*selectionParams.nGM)); % update waitbar
    end
    
    % check whether results are within tolerance, and stop optimization if so
    if within_tolerance(sampleSmall, targetSa, selectionParams)
        break;
    end    
end

close(hw); % close waitbar

% Save final selection for output
IMs.sampleSmall = sampleSmall;


delete(parobj);
end

function [scaleFac, minDev] = bestScaleFactorPar(sampleBig,sampleSmall,meanReq,sigma,weights,maxScale)
% Identifies the best scaled ground motions to be used with the greedy
% algortihm

% Determine size of sampleSmall for standard deviation calculations
scales = 0.1:0.1:maxScale;
scaleFac = zeros(size(sampleBig,1),1);
[nGM,~] = size(sampleSmall);
newGM = nGM+1;
minDev = zeros(length(scales),1);

parfor i = 1:size(sampleBig,1)
    devTotal = zeros(length(scales),1);
    
    for j=1:length(scales)
        
        sampleSmallNew = [sampleSmall;sampleBig(i,:)+log(scales(j))];
        
        % Compute deviations from target
        avg = sum(sampleSmallNew)./newGM;
        devMean = avg - meanReq;
        devSig = sqrt((1/(nGM))*sum((sampleSmallNew-repmat(avg,newGM,1)).^2))-sigma;
        devTotal(j) = weights(1) * sum(devMean.^2) + weights(2) * sum(devSig.^2);
        
    end
    
    [minDev(i), minID] = min(devTotal);
    
    if minDev(i) == 100000;
        scaleFac(i) = -99;
    else
        scaleFac(i) = scales(minID);
    end
    
    
end
end



function [ devTotal ] = ParLoop( devTotal, scaleFac, selectionParams, sampleSmall, IMs, meanReq, stdevs, emp_cdf )
% Parallel loop to use within greedy optimization
optType = selectionParams.optType;
if all(devTotal) && optType == 0 
    return;
end

TgtPer = selectionParams.TgtPer;
cond = selectionParams.cond;
isScaled = selectionParams.isScaled;
weights = selectionParams.weights;
penalty = selectionParams.penalty;
maxScale = selectionParams.maxScale;
recID = IMs.recID;
sampleBig = IMs.sampleBig;



parfor j = 1:selectionParams.nBig
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
        [devTotal(j)] = KS_stat(TgtPer, emp_cdf, sampleSmallTemp, meanReq, stdevs);
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


