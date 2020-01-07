function [ IMs ] = find_ground_motionsV( selectionParams, simulatedSpectra, IMs )
% Find best matches to the simulated spectra from ground-motion database

% Simplify notation
nPer = length(selectionParams.TgtPer);
nPerV = length(selectionParams.TgtPerV);

% Extract data for V component of GM before resetting variable for H
% component
simulatedSpectraV = simulatedSpectra(:,(nPer+1):end);
simulatedSpectra = simulatedSpectra(:,1:nPer);

% Relative importance of components of GM
wH = 1-selectionParams.weightV;

% Identify periods for scaling
if selectionParams.cond == 1 
    scaleFacIndex = selectionParams.indTcond;
else
    scaleFacIndex = (1:nPer)';
end
scaleFacIndexV = (1:nPerV)'; % If sepScaleV=1, V components scaled to target at selectionParams.TgtPerV

% Initialize vectors
IMs.recID = zeros(selectionParams.nGM,1);
IMs.sampleSmall = [];
IMs.scaleFac = ones(selectionParams.nGM,1);
IMs.sampleSmallV = [];
IMs.scaleFacV = ones(selectionParams.nGM,1);

% Find database spectra most similar to each simulated spectrum
for i = 1:selectionParams.nGM % for each simulated spectrum
    err = zeros(selectionParams.nBig,1); % initialize error matrix
    scaleFac = ones(selectionParams.nBig,1); % initialize scale factors to 1 (these are never changed if no scaling is allowed)
    scaleFacV = ones(selectionParams.nBig,1); % initialize scale factors to 1 (these are never changed if no scaling is allowed)
    
    % compute scale factors and errors for each candidate ground motion
    for j=1:selectionParams.nBig
        if selectionParams.isScaled % if scaling is allowed, compute scale factor (otherwise it is already defined as equal to 1)            
            scaleFac(j) = geomean( simulatedSpectra(i,scaleFacIndex) ./ exp(IMs.sampleBig(j,scaleFacIndex)) );
            if selectionParams.sepScaleV == 1 % if separate scale factors for V component desired, compute based on target spectrum for V component
                scaleFacV(j) = geomean( simulatedSpectraV(i,scaleFacIndexV) ./ exp(IMs.sampleBigV(j,scaleFacIndexV)) );
            else
                scaleFacV(j) = scaleFac(j);
            end
        end
        err(j) = sum((log(exp(IMs.sampleBig(j,:))*scaleFac(j)) - log(simulatedSpectra(i,:))).^2)*wH + sum((log(exp(IMs.sampleBigV(j,:))*scaleFacV(j)) - log(simulatedSpectraV(i,:))).^2)*(1-wH);
    end
    err(IMs.recID(1:i-1)) = 1e6; % exclude previously-selected ground motions 
    err(scaleFac > selectionParams.maxScale | scaleFac < (1/selectionParams.maxScale)) = 1e6; % exclude ground motions requiring too large of a scale factor
    err(scaleFacV > selectionParams.maxScale | scaleFacV < (1/selectionParams.maxScale)) = 1e6; % exclude ground motions requiring too large of a scale factor
   
    % find minimum-error ground motion    
    [tmp, IMs.recID(i)] = min(err);
    assert(tmp < 1e6, 'Warning: problem with simulated spectrum. No good matches found');
    IMs.scaleFac(i) = scaleFac(IMs.recID(i)); % store scale factor
    IMs.sampleSmall = [IMs.sampleSmall; log(exp(IMs.sampleBig(IMs.recID(i),:))*scaleFac(IMs.recID(i)))]; % store scaled log spectrum
    IMs.scaleFacV(i) = scaleFacV(IMs.recID(i)); % store scale factor
    IMs.sampleSmallV = [IMs.sampleSmallV; log(exp(IMs.sampleBigV(IMs.recID(i),:))*scaleFacV(IMs.recID(i)))]; % store scaled log spectrum
    
end


end

