function [ scaleIndex, corrMatrix, Targets, optInputs ] = ComputeTargets( indPer, knownPer, sa, sigma, useVar, eps_bar, optInputs )
% ComputeTargets will calculate and return the target mean spectrum, target
% covariance matrix, and target correlation matrix for ground motion
% selection. The index/indicies of PerTgt that will need to be scaled is
% also returned. The predicted spectral accelerations and their
% corresponding standard deviations can be computed using any GMPE the user
% chooses. The function accepts the following arguments:
%           indPer          : periods at which target means will be
%                             calculated
%           knownPer        : available periods from the database
%           sa              : spectral accelerations from the GMPE
%           sigma           : standard deviations from the GMPE
%           useVar          : user input for calculating variance
%           eps_bar         : user input for epsilon (conditional
%                             selection)
%           optInputs       : optimization user inputs 
%% Compute target mean spectrum from results of Campbell and Bozorgnia GMPE
% Initialize variables for computing targets 
% perKnownRec = find(perKnownCorr == optInputs.T1);

% (Log) Response Spectrum Mean: meanReq
% Define indicies at which spectra will be scaled
if optInputs.cond == 1 
    scaleIndex = optInputs.indT1; 
    
%     % remove values at T1 in order to compute mean values to match known
%     % periods
%     if ~any(knownPer == optInputs.T1)
%         sa(perKnownRec) = [];
%         sigma(perKnownRec) = [];
%     end
    
    % compute correlations and the conditional mean spectrum
    rho = zeros(1,length(indPer));  
    for i = 1:length(indPer)
        rho(i) = baker_jayaram_correlation(knownPer(indPer(i)), optInputs.TgtPer(optInputs.indT1));
    end
    Targets.meanReq = log(sa(indPer)) + sigma(indPer).*eps_bar.*rho;
        
    % define the spectral accleration at T1 that all ground motions will be
    % scaled to
    % add to main script
%     optInputs.lnSa1 = Targets.meanReq(optInputs.indT1);
    
elseif optInputs.cond == 0
    scaleIndex = (1:length(optInputs.TgtPer))';
    Targets.meanReq = log(sa(indPer));
end


%% Compute covariances at all periods 
covMatrix = zeros(length(sa));
corrMatrix = zeros(length(sa));
for i=1:length(sa)
    for j=1:length(sa)
        % Periods
        Ti = knownPer(i);
        Tj = knownPer(j);

        % Means and variances
        varT = sigma(optInputs.indT1)^2;
        sigma22 = varT;
        var1 = sigma(i)^2;
        var2 = sigma(j)^2;

        if optInputs.cond == 1 
            sigma11 = [var1 baker_jayaram_correlation(Ti, Tj)*sqrt(var1*var2);baker_jayaram_correlation(Ti, Tj)*sqrt(var1*var2) var2];
            sigma12 = [baker_jayaram_correlation(Ti, optInputs.T1)*sqrt(var1*varT);baker_jayaram_correlation(optInputs.T1, Tj)*sqrt(var2*varT)];
            sigmaCond = sigma11 - sigma12*inv(sigma22)*(sigma12)';
            % Correlations
            covMatrix(i,j) = sqrt(sigmaCond(1,1)*sigmaCond(2,2));
            corrMatrix(i,j) = sigmaCond(1,2)/covMatrix(i,j);
        elseif optInputs.cond == 0
            % Correlations
            covMatrix(i,j) = sqrt(var1*var2);
            corrMatrix(i,j) = baker_jayaram_correlation(Ti, Tj);
        end

    end
end

%% Compute target conditional coveriance at periods of interest 
if useVar == 0
    Targets.covReq = zeros(length(optInputs.TgtPer));
else
    Targets.covReq = corrMatrix(indPer,indPer).*covMatrix(indPer,indPer);
    
    % for conditional selection only, ensure that variances will be zero at
    % all values of T1 (but not exactly 0.0, for MATLAB spectra simulations)
    if optInputs.cond == 1
        Targets.covReq(optInputs.indT1,:) = repmat(1e-17,1,length(optInputs.TgtPer));
        Targets.covReq(:,optInputs.indT1) = repmat(1e-17,length(optInputs.TgtPer),1);
    end    

end

