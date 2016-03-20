function [ IMs ] = find_ground_motions( optInputs, simulatedSpectra, IMs )
% Find best matches to the simulated spectra from ground-motion database


if optInputs.cond == 1 
    scaleFacIndex = optInputs.indT1;
else
    scaleFacIndex = (1:length(optInputs.TgtPer))';
end

% Initialize vectors
IMs.recID = zeros(optInputs.nGM,1);
IMs.sampleSmall = [];
IMs.scaleFac = ones(optInputs.nGM,1);

% Find database spectra most similar to each simulated spectrum
for i = 1:optInputs.nGM % for each simulated spectrum
    err = zeros(optInputs.nBig,1); % initialize error matrix
    scaleFac = ones(optInputs.nBig,1); % initialize scale factors to 1 (these are never changed if no scaling is allowed)
    
    % compute scale factors and errors for each candidate ground motion
    for j=1:optInputs.nBig
        if optInputs.isScaled % if scaling is allowed, compute scale factor (otherwise it is already defined as equal to 1)
            scaleFac(j) = sum(exp(IMs.sampleBig(j,scaleFacIndex)).*simulatedSpectra(i,scaleFacIndex))/sum(exp(IMs.sampleBig(j,scaleFacIndex)).^2);
        end
        err(j) = sum((log(exp(IMs.sampleBig(j,:))*scaleFac(j)) - log(simulatedSpectra(i,:))).^2);
    end
    err(IMs.recID(1:i-1)) = 1000000; % exclude previously-selected ground motions 
    err(scaleFac > optInputs.maxScale) = 1000000; % exclude ground motions requiring too large of a scale factor
   
    % find minimum-error ground motion    
    [tmp, IMs.recID(i)] = min(err);
    assert(tmp < 1000000, 'Warning: problem with simulated spectrum. No good matches found');
    IMs.scaleFac(i) = scaleFac(IMs.recID(i)); % store scale factor
    IMs.sampleSmall = [IMs.sampleSmall; log(exp(IMs.sampleBig(IMs.recID(i),:))*scaleFac(IMs.recID(i)))]; % store scaled log spectrum
end


end

