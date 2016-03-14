% Identifies the best scaled ground motions to be used with the greedy
% algortihm

function [scaleFac, minDev] = bestScaleFactor(sampleBig,sampleSmall,meanReq,sigma,weights,maxScale)

% Determine size of sampleSmall for standard deviation calculations
scales = 0.1:0.1:maxScale;
scaleFac = zeros(size(sampleBig,1),1);
[nGM,~] = size(sampleSmall);
newGM = nGM+1;
devTotal = zeros(length(scales),1);
minDev = zeros(length(scales),1);
% minDevFinal = zeros(size(sampleBig,1),1);

for i = 1:size(sampleBig,1)
    for j=1:length(scales)
        
        sampleSmallNew = [sampleSmall;sampleBig(i,:)+log(scales(j))];
        
        % Compute deviations from target
        avg = sum(sampleSmallNew)./newGM;
        devMean = avg - meanReq;
        devSig = sqrt((1/(newGM-1))*sum((sampleSmallNew-repmat(avg,newGM,1)).^2))-sigma;
        devTotal(j) = weights(1) * sum(devMean.^2) + weights(2) * sum(devSig.^2);
        
    end
    
    [minDev(i), minID] = min(devTotal);
    
    if minDev(i) == 100000;
        scaleFac(i) = -99;
    else
        scaleFac(i) = scales(minID);
    end

end

