function [ TgtMean, TgtCovs ] = ComputeTargets( RotD, arb, indPer, knownPer, useVar, eps_bar, optInputs, M_bar, Rjb, Fault_Type, region, z1, Vs30)
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

% compute the median and standard deviations of RotD50 response spectrum values 
% if the GMPE arguments are not available, input the predicted sa and sigma
% values 
if nargin < 8
    sa = []; % input median results from any other GMPE (of length knownPer <= 10)
    sigma = []; % input predicted sigmas
else
    [sa, sigma] = BSSA_2014_nga(M_bar, knownPer(knownPer<=10), Rjb, Fault_Type, region, z1, Vs30);
end

% modify spectral targets if RotD100 values were specified for
% two-component selection, see the following document for more details:
%  Shahi, S. K., and Baker, J. W. (2014). "NGA-West2 models for ground-
%  motion directionality." Earthquake Spectra, 30(3), 1285-1300.
if RotD == 100 && arb == 2 
   [ rotD100Ratio, rotD100Sigma ] = SB_2014_ratios( knownPer ); 
   sa = sa.*rotD100Ratio; % from equation (1) of the above paper
   sigma = sqrt ( sigma.^2 + rotD100Sigma .^2); 
end

% (Log) Response Spectrum Mean: TgtMean
if optInputs.cond == 1   
    % compute correlations and the conditional mean spectrum
    rho = zeros(1,length(sa));  
    for i = 1:length(sa)
        rho(i) = baker_jayaram_correlation(knownPer(i), optInputs.TgtPer(optInputs.indT1));
    end
    TgtMean = log(sa) + sigma.*eps_bar.*rho;

elseif optInputs.cond == 0
    TgtMean = log(sa);
end


%% Compute covariances and correlations at all periods 
TgtCovs = zeros(length(sa));
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
            sigmaCorr = baker_jayaram_correlation(Ti,Tj)*sqrt(var1*var2);
            sigma11 = [var1 sigmaCorr;sigmaCorr var2];
            sigma12 = [baker_jayaram_correlation(Ti, optInputs.T1)*sqrt(var1*varT);baker_jayaram_correlation(optInputs.T1, Tj)*sqrt(var2*varT)];
            sigmaCond = sigma11 - sigma12*inv(sigma22)*(sigma12)';
            
            TgtCovs(i,j) = sigmaCond(1,2);
        elseif optInputs.cond == 0
            TgtCovs(i,j) = baker_jayaram_correlation(Ti, Tj)*sqrt(var1*var2);
        end

    end
end

%% Compute target conditional coveriance at periods of interest 
if useVar == 0
    TgtCovs = zeros(length(knownPer));
else
    % for conditional selection only, ensure that variances will be zero at
    % all values of T1 (but not exactly 0.0, for MATLAB spectra simulations)
    if optInputs.cond == 1
        TgtCovs(indPer(optInputs.indT1),:) = 1e-17;
        TgtCovs(:,indPer(optInputs.indT1)) = 1e-17;
    end
end

