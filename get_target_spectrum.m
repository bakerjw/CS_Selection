function [ Tgts ] = get_target_spectrum( RotD, arb, knownPer, useVar, eps_bar, optInputs, indPer, rup)
% Calculate and return the target mean spectrum, target and covariance
% matrix at available periods. Predicted spectral accelerations and
% corresponding standard deviations are computed using gmpeBSSA_2014 but
% can be replaced by a user-defined GMPE. If the GMPE is altered then the
% input arguments will likely also need to be adjusted
%
% INPUTS
%           RotD            : =50 or =100
%           arb             : =1 for single-component or =2 for two-component
%           knownPer        : available periods from the database
%           indPer          : period indicies of target response spectrum
%           useVar          : user input for calculating variance
%           eps_bar         : user input for epsilon (conditional
%                             selection)
%           optInputs       : optimization user inputs 
%           rup             :  A structure with parameters that specify the rupture scenario
%                              for the purpose of evaluating a GMPE
%
%               .M_bar           : earthquake magnitude
%               .Rjb             : closest distance to surface projection of the fault rupture (km)
%               .Fault_Type      : =0 for unspecified fault
%                                   =1 for strike-slip fault
%                                   =2 for normal fault
%                                   =3 for reverse fault
%               .region          : =0 for global (incl. Taiwan)
%                                   =1 for California 
%                                   =2 for Japan 
%                                   =3 for China or Turkey 
%                                   =4 for Italy
%               .z1              : basin depth (km); depth from ground surface to the 1km/s shear-wave horizon, =999 if unknown
%               .Vs30            : average shear wave velocity in the top 30m of the soil (m/s)
%
%
% Outputs (these could be replaced by user-specified matrices if desired
%           TgtMean         : (log) target mean response spectrum with length of knownPer
%           TgtCovs         : target covariance matrix with size [length(knownPer) length(knownPer)]


%% Compute target mean spectrum
% compute the median and standard deviations of RotD50 response spectrum values 
[sa, sigma] = gmpe_bssa_2014(rup.M_bar, knownPer(knownPer<=10), rup.Rjb, rup.Fault_Type, rup.region, rup.z1, rup.Vs30);

% modify spectral targets if RotD100 values were specified for
% two-component selection, see the following document for more details:
%  Shahi, S. K., and Baker, J. W. (2014). "NGA-West2 models for ground-
%  motion directionality." Earthquake Spectra, 30(3), 1285-1300.
if RotD == 100 && arb == 2 
   [ rotD100Ratio, rotD100Sigma ] = gmpe_sb_2014_ratios( knownPer ); 
   sa = sa.*rotD100Ratio; % from equation (1) of the above paper
   sigma = sqrt ( sigma.^2 + rotD100Sigma .^2); 
end

% back-calculate and epsilon to get the targe Sa(T_cond) if needed
if 1==1
%     eps_bar = (log(Sa_star) - log(Sa_1))/sigma_1;
end

% (Log) Response Spectrum Mean: TgtMean
if optInputs.cond == 1   
    % compute correlations and the conditional mean spectrum
    rho = zeros(1,length(sa));  
    for i = 1:length(sa)
        rho(i) = gmpe_bj_2008_corr(knownPer(i), optInputs.TgtPer(optInputs.indT1));
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
            sigma12 = [gmpe_bj_2008_corr(Ti, optInputs.T1)*sqrt(var1*varT);baker_jayaram_correlation(optInputs.T1, Tj)*sqrt(var2*varT)];
            sigmaCond = sigma11 - sigma12*inv(sigma22)*(sigma12)';
            TgtCovs(i,j) = sigmaCond(1,2);
        elseif optInputs.cond == 0
            TgtCovs(i,j) = baker_jayaram_correlation(Ti, Tj)*sqrt(var1*var2);
        end
    end
end

% over-write coveriance matrix with zeros if no variance is desired in the ground motion selection
if useVar == 0
    TgtCovs = zeros(size(TgtCovs));
end

% Store target mean and covariance matrix at target periods
Tgts.meanReq = TgtMean(indPer); 
Tgts.covReq  = TgtCovs(indPer,indPer);
Tgts.stdevs  = sqrt(diag(Tgts.covReq))';


