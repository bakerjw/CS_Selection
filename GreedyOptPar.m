function [ sampleSmall, finalRecords, finalScaleFactors ] = GreedyOptPar( optInputs, Tgts, IMs )
% This function will perform a greedy optimization on a set of ground
% motions using the sum of squared errors approach to check the set of
% selected ground motions against target means and variances using
% parallelization from MATLAB's "parfor" function

% Redefine relevant variables from structures
% IMs structure
sampleSmall = IMs.sampleSmall;
sampleBig   = IMs.sampleBig;

% optInputs structure
isScaled    = optInputs.isScaled;
cond        = optInputs.cond;
nGM         = optInputs.nGM;
penalty     = optInputs.penalty;
weights     = optInputs.weights;
maxScale    = optInputs.maxScale;
nBig        = optInputs.nBig;
recID       = optInputs.recID;
rec         = optInputs.rec;

% Tgts structure
meanReq     = Tgts.meanReq;
stdevs      = Tgts.stdevs;

% Define scale factors for conditional selection
if cond == 1
    lnSa1       = optInputs.lnSa1;
    sampleBig_scale = sampleBig(:,rec);
    condScale = exp(lnSa1)./exp(sampleBig_scale);
end

% Define new sampleSmall matrices and recID vectors that take away each
% respective ground motion spectrum at a time
sampleSmall_nGM = zeros(nGM-1,size(sampleSmall,2),nGM);
recID_nGM = zeros(nGM-1,nGM);
for i = 1:nGM
    sampleSmall_nGM(:,:,i) = sampleSmall([1:i-1, i+1:end],:);
    recID_nGM(:,i) = recID([1:i-1, i+1:end]);
end

% open parallel pool workers
numWorkers = 2;
parobj = parpool(numWorkers);

% Initialize final variables (not temporary, check for temporary variables
% that will disappear after parfor)
sampleSmall_2 = zeros(nGM,size(sampleSmall,2),nGM);
recID_2 = zeros(nGM,nGM);
finalScaleFac = zeros(nGM);
% parallelize optimization for multiple ground motions
devTotal = zeros(nBig,1);
sampleSmall_1 = sampleSmall_nGM(:,:,i);
recID_1 = recID_nGM(:,i);

if isScaled == 1
    if cond == 0
        [scaleFac, devTotal] = bestScaleFactor(sampleBig, sampleSmall_1, meanReq, stdevs, weights, maxScale);
    elseif cond == 1
        scaleFac = condScale;
    end
end

% Try to add a new spectra to the subset list
for j = 1:nBig
    if scaleFac(j) > maxScale
        devTotal(j) = 100000;
    else
        sampleSmallTemp = [sampleSmall_1;sampleBig(j,:)+log(scaleFac(j))];
        
        % Calculate the appropriate measure of deviation and store in
        % devTotal
        if cond == 1 || (cond == 0 && isScaled == 0)
            % Compute deviations from target
            devMean = mean(sampleSmallTemp) - meanReq;
            devSig = std(sampleSmallTemp) - stdevs;
            devTotal(j) = weights(1) * sum(devMean.^2) + weights(2) * sum(devSig.^2);
        end
        % Penalize bad spectra (set penalty to zero if this is not required)
        if penalty ~= 0
            for m=1:size(sampleSmallTemp,1)
                devTotal(j) = devTotal(j) + sum(abs(exp(sampleSmallTemp(m,:))>exp(meanReq+3*stdevs'))) * penalty;
            end
        end
        
        % Should cause improvement and record should not be repeated
        if (any(recID_1 == j))
            devTotal(j) = 100000;
        end
    end
end
[~, minID] = min(devTotal);
% Add new element in the right slot
if isScaled == 1
    finalScaleFac(i) = scaleFac(minID);
else
    finalScaleFac(i) = 1;
end
sampleSmallNew = sampleSmall_nGM(:,:,i);
recIDNew = recID_nGM(:,i);
sampleSmallFinal = [sampleSmallNew(1:i-1,:);sampleBig(minID,:)+log(scaleFac(minID));sampleSmallNew(i:end,:)];
recIDFinal = [recIDNew(1:i-1);minID;recIDNew(i:end)];

sampleSmall_2(:,:,i) = sampleSmallFinal;
recID_2(:,i) = recIDFinal;



% Output information
sampleSmall = sampleSmall_2;
finalRecords = recID_2;
finalScaleFactors = finalScaleFac';


delete(parobj);

end

