function [scaleFac, minDev] = compute_scale_factor(sampleBig,sampleSmall,meanReq,sigma,weights,maxScale)
% Compute scale factors for a set of candidate ground motions, when using 
% an Unconditional Spectrum target 


nGM = size(sampleSmall,1) + 1; % number of ground motions
scales = 0.1:0.1:maxScale;
scaleFac = zeros(size(sampleBig,1),1);
devTotal = zeros(length(scales),1);
minDev = zeros(length(scales),1);

for i = 1:size(sampleBig,1)
    for j=1:length(scales)
        
        sampleSmallNew = [sampleSmall;sampleBig(i,:)+log(scales(j))]; % candidate with ground motion i and scale factor j
        
        % Compute deviations from target
        avg = sum(sampleSmallNew)./nGM;
        devMean = avg - meanReq;
        devSig = sqrt((1/(nGM-1))*sum((sampleSmallNew-repmat(avg,nGM,1)).^2))-sigma;
        devTotal(j) = weights(1) * sum(devMean.^2) + weights(2) * sum(devSig.^2);
        
    end
    
    [minDev(i), minID] = min(devTotal);
    scaleFac(i) = scales(minID);

end

