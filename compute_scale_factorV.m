function [scaleFac, minDev] = compute_scale_factorV(IMs,sampleSmall,sampleSmallV,targetSa,selectionParams)
% Compute scale factors for a set of candidate ground motions, when using 
% an Unconditional Spectrum target 

% Extract input data from structures
nPer = length(selectionParams.TgtPer);
meanReq = targetSa.meanReq(1,1:nPer);
meanReqV = targetSa.meanReq(1,(nPer+1):end);
sigma = targetSa.stdevs(1,1:nPer);
sigmaV = targetSa.stdevs(1,(nPer+1):end);
sampleBig = IMs.sampleBig;
sampleBigV = IMs.sampleBigV;
weights = selectionParams.weights;
maxScale = selectionParams.maxScale;
wH = 1-selectionParams.weightV;

nGM = size(sampleSmall,1) + 1; % number of ground motions
scales = (1/maxScale):0.1:maxScale;
scaleFac = zeros(size(sampleBig,1),1);
devTotal = zeros(length(scales),1);
minDev = zeros(length(scales),1);

for i = 1:size(sampleBig,1)
    for j=1:length(scales)
        
        sampleSmallNew = [sampleSmall;sampleBig(i,:)+log(scales(j))]; % candidate with ground motion i and scale factor j
        sampleSmallNewV = [sampleSmallV;sampleBigV(i,:)+log(scales(j))]; % candidate with ground motion i and scale factor j
        
        % Compute deviations from target
        avg = sum(sampleSmallNew)./nGM;
        avgV = sum(sampleSmallNewV)./nGM;
        devMean = avg - meanReq;   
        devMeanV = avgV - meanReqV;      
        devSig = sqrt((1/(nGM-1))*sum((sampleSmallNew-repmat(avg,nGM,1)).^2))-sigma;
        devSigV = sqrt((1/(nGM-1))*sum((sampleSmallNewV-repmat(avgV,nGM,1)).^2))-sigmaV;
        devTotal(j) = weights(1) * (sum(devMean.^2)*wH + sum(devMeanV.^2)*(1-wH)) + weights(2) * (sum(devSig.^2)*wH + sum(devSigV.^2)*(1-wH));
        
    end
    
    [minDev(i), minID] = min(devTotal);
    scaleFac(i) = scales(minID);

end

