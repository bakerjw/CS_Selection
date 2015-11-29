function [ scaleIndex, corrMatrix, Targets ] = ComputeTargets( recPer, perKnown, perKnownCorr, ... 
                                                            saCorr, sigmaCorr, useVar, eps_bar, optInputs )
%%

perKnownRec = find(perKnownCorr == optInputs.T1);
sa = saCorr;
sigma = sigmaCorr; 

% Estimate target means and covariances
% (Log) Response Spectrum Mean: meanReq
if optInputs.cond == 1 
    % define the index needed for scaling
    scaleIndex = optInputs.rec; % optInputs.T1index
    
    % remove values at T1 in order to compute mean values to match known
    % periods
    if ~any(perKnown == optInputs.T1)
    sa(perKnownRec) = []; 
    sigma(perKnownRec) = [];
    end
    
    % compute correlations and the conditional mean spectrum
    rho = zeros(1,length(recPer));  
    for i = 1:length(recPer)
        rho(i) = baker_jayaram_correlation(perKnown(recPer(i)), optInputs.T1);
    end
    Targets.meanReq = log(sa(recPer)) + sigma(recPer).*eps_bar.*rho;
        
    % define the spectral accleration at T1 that all ground motions will be
    % scaled to
    optInputs.lnSa1 = Targets.meanReq(optInputs.rec);
    
elseif optInputs.cond == 0
    scaleIndex = (1:length(optInputs.PerTgt))';
    Targets.meanReq = log(sa(recPer));
end

%% Estimate covariances at all available periods from the Baker and Jayaram (2008) model
covReqPart = zeros(length(perKnownCorr));
corrMatrix = zeros(length(perKnownCorr));
for i=1:length(perKnownCorr)
    for j=1:length(perKnownCorr)
        % Periods
        Ti = perKnownCorr(i);
        Tj = perKnownCorr(j);

        % Means and variances
        varT = sigmaCorr(optInputs.rec)^2;
        sigma22 = varT;
        var1 = sigmaCorr(i)^2;
        var2 = sigmaCorr(j)^2;

        if optInputs.cond == 1 
            sigma11 = [var1 baker_jayaram_correlation(Ti, Tj)*sqrt(var1*var2);baker_jayaram_correlation(Ti, Tj)*sqrt(var1*var2) var2];
            sigma12 = [baker_jayaram_correlation(Ti, optInputs.T1)*sqrt(var1*varT);baker_jayaram_correlation(optInputs.T1, Tj)*sqrt(var2*varT)];
            sigmaCond = sigma11 - sigma12*inv(sigma22)*(sigma12)';
            % Correlations
            covReqPart(i,j) = sqrt(sigmaCond(1,1)*sigmaCond(2,2));
            corrMatrix(i,j) = sigmaCond(1,2)/sqrt(sigmaCond(1,1)*sigmaCond(2,2));
        elseif optInputs.cond == 0
            % Correlations
            covReqPart(i,j) = sqrt(var1*var2);
            corrMatrix(i,j) = baker_jayaram_correlation(Ti, Tj);
        end

    end
end

%% (Log) Response Spectrum Covariance: covReq
if useVar == 0
    Targets.covReq = zeros(length(optInputs.PerTgt));
else
    if optInputs.cond == 1 && ~any(perKnown == optInputs.T1)
        corrReqNew = corrMatrix;
        corrReqNew(perKnownRec,:) = []; corrReqNew(:,perKnownRec) = [];
        covReqPartNew = covReqPart;
        covReqPartNew(perKnownRec,:) = []; covReqPartNew(:,perKnownRec) = [];
        Targets.covReq = corrReqNew(recPer,recPer).*covReqPartNew(recPer,recPer);
    else
        Targets.covReq = corrMatrix(recPer,recPer).*covReqPart(recPer,recPer);
    end
    
    % for conditional selection only, ensure that variances will be zero at
    % all values of T1 (but not exactly 0.0, for MATLAB spectra simulations)
    if optInputs.cond == 1
        Targets.covReq(optInputs.rec,:) = repmat(1e-17,1,length(optInputs.PerTgt));
        Targets.covReq(:,optInputs.rec) = repmat(1e-17,length(optInputs.PerTgt),1);
    end    

end

