function [ TgtMean, TgtCovs ] = ComputeTargets( RotD, arb, indPer, knownPer, useVar, eps_bar, optInputs, M_bar, Rjb, Fault_Type, region, z1, Vs30)
%% Calculate and return the target mean spectrum, target and covariance
% matrix at available periods. Predicted spectral accelerations and
% corresponding standard deviations are computed using BSSA_2014_nga but
% can be replaced by a user-defined GMPE. 
%           RotD            : =50 or =100
%           arb             : =1 for single-component or =2 for two-component
%           knownPer        : available periods from the database
%           indPer          : period indicies of target response spectrum
%           useVar          : user input for calculating variance
%           eps_bar         : user input for epsilon (conditional
%                             selection)
%           optInputs       : optimization user inputs 
%           M_bar           : earthquake magnitude
%           Rjb             : closest distance to surface projection of the fault rupture (km)
%           Fault_Type      : =0 for unspecified fault
%                             =1 for strike-slip fault
%                             =2 for normal fault
%                             =3 for reverse fault
%           region          : =0 for global (incl. Taiwan)
%                             =1 for California 
%                             =2 for Japan 
%                             =3 for China or Turkey 
%                             =4 for Italy
%           z1              : basin depth (km); depth from ground surface to the 1km/s shear-wave horizon, =999 if unknown
%           Vs30            : average shear wave velocity in the top 30m of the soil (m/s)
%% Compute target mean spectrum
% compute the median and standard deviations of RotD50 response spectrum values 
[sa, sigma] = BSSA_2014_nga(M_bar, knownPer(knownPer<=10), Rjb, Fault_Type, region, z1, Vs30);

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
        
        % Covariances 
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
    TgtCovs = zeros(length(sa));
else
    % for conditional selection only, ensure that variances will be zero at
    % all values of T1 (but not exactly 0.0, for MATLAB spectra simulations)
    if optInputs.cond == 1
        TgtCovs(indPer(optInputs.indT1),:) = 1e-17;
        TgtCovs(:,indPer(optInputs.indT1)) = 1e-17;
    end
end

