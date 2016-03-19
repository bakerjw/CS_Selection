% This code is used to select ground motions with response spectra 
% representative of a target scenario earthquake, as predicted by a ground 
% motion model. Spectra can be selected to be consistent with the full
% distribution of response spectra, or conditional on a spectral amplitude
% at a given period (i.e., using the conditional spectrum approach). 
% Single-component or two-component motions can be selected, and several
% ground motion databases are provided to search in. Further details are
% provided in the following documents:
%
%   N. Jayaram, T. Lin and and Baker, J. W. (2011). A computationally
%   efficient ground-motion selection algorithm for matching a target
%   response spectrum mean and variance, Earthquake Spectra, 27(3), 797-815.
%
% created by Nirmal Jayaram, Ting Lin, Jack W. Baker, Official release 7 June, 2010 
%
% modified by Cynthia Lee and Jack Baker, last updated, 14 March, 2016
%% OUTPUT VARIABLES
%
% finalRecords      : Record numbers of selected records
% finalScaleFactors : Corresponding scale factors
%
% (these variables are also output to a text file specified by the outputFile variable)
%
% The final cell in this m file shows how to plot the selected spectra
% using this information. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variable definitions and initial set of user inputs
% Specify a database with a needed ground motion data. Provided databases
% include 'NGA_W1_meta_data.mat', 'NGA_W2_meta_data.mat',
% 'BBP_GP_meta_data.mat', 'BBP_SDSU_meta_data.mat',
% and'BBP_EXSIM_meta_data.mat' Further documentation of these databases can
% be found at 'WorkspaceDocumentation***.txt'.
%
% Variable definitions for loading data:
% cond         : 0 to run unconditional selection
%                1 to run conditional
% arb          : 1 for single-component selection and arbitrary component sigma
%                2 for two-component selection and average component sigma
% RotD         : 50 to use SaRotD50 data
%              : 100 to use SaRotD100 data
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Certain variables are stored in data structures which are later used in
% the optimization function. Required user input values are indicated for
% the user. Some variables defined here are calculated within this script
% or other functions. The data structures are as follows:
%
% optInputs: input values needed to run the optimization function
%            isScaled   : The user will input 1 to allow records to be 
%                         scaled, and input 0 otherwise 
%            maxScale   : The maximum allowable scale factor
%            tol        : User input percent error tolerance to determine
%                         whether or not optimization can be skipped (only
%                         used for SSE optimization)
%            optType    : For greedy optimization, the user will input a 0
%                         to use the sum of squared errors approach to 
%                         optimize the selected spectra, or a 1 to use 
%                         D-statistic calculations from the KS-test
%            penalty    : If a penalty needs to be applied to avoid selecting
%                         spectra that have spectral acceleration values 
%                         beyond 3 sigma at any of the periods, set a value
%                         here. Use 0 otherwise.
%            weights    : [Weights for error in mean, standard deviation 
%                         and skewness] e.g., [1.0,2.0 0.3] 
%            nLoop      : Number of loops of optimization to perform.
%                         Default value = 2
%            nBig       : The number of spectra that will be searched
%            indT1      : This is the index of T1, the conditioning period
%            recID      : This is a vector of index values for chosen
%                         spectra
% 
% Tgts    :  The target values (means and covariances) being matched
%            meanReq    : Estimated target response spectrum means (vector of
%                         logarithmic spectral values, one at each period)
%            covReq     : Matrix of response spectrum covariances
%            stdevs     : A vector of standard deviations at each period
% 
% IMs     :  The intensity measure values (from SaKnown) chosen and the 
%            values available
%            sampleSmall: matrix of selected logarithmic response spectra 
%            sampleBig  : The matrix of logarithmic spectra that will be 
%                          searched
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If a database other than the provided databases is used, define the
% following variables:
% SaKnown   : (N*P matrix)
%             This is a matrix of Sa values at different periods (P) for
%             available ground-motion time histories (N).
% knownPer  : The set of P known periods.
% 
% soil_Vs30        : Soil Vs30 values corresponding to all the records
% magnitude        : Magnitude of all the records
% closest_D        : Closest distance for all the records
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User inputs begin here
% Ground motion database and type of selection 
databaseFile         = 'NGA_W2_meta_data'; % filename of the target database
optInputs.cond       = 1;
arb                  = 2; 
RotD                 = 50; 

% Number of ground motions and spectral periods of interest
optInputs.nGM        = 10;      % number of ground motions to be selected 
optInputs.T1         = 0.5;     % Period at which spectra should be scaled and matched 
Tmin                 = 0.1;     % smallest spectral period of interest
Tmax                 = 10;      % largest spectral period of interest
optInputs.TgtPer     = logspace(log10(Tmin),log10(Tmax),30); % Periods at which the target spectrum needs to be computed (logarithmically spaced)

% other parameters to scale motions and evaluate selections 
optInputs.isScaled   = 1;       
optInputs.maxScale   = 4;       
optInputs.tol        = 0; 
optInputs.optType    = 1; 
optInputs.penalty    = 0;
optInputs.weights    = [1.0 2.0 0.3];
optInputs.nLoop      = 2;

% User inputs to specify the target earthquake scenario
M_bar       = 7.52;     % earthquake magnitude
R_bar       = 11.6;     % distance corresponding to the target scenario earthquake
Rjb         = R_bar;    % closest distance to surface projection of the fault rupture (km)
eps_bar     = 1.9524;   % epsilon value (used only for conditional selection)
                        % BackCalcEpsilon is a supplementary script that
                        % will back-calculate epsilon values based on M, R,
                        % and Sa inputs from USGS deaggregations
Vs30        = 259;      % average shear wave velocity in the top 30m of the soil (m/s)
z1          = 999;      % basin depth (km); depth from ground surface to the 1km/s shear-wave horizon,
                        % =999 if unknown
useVar      = 1;        % =1 to use ground motion model variance, =0 to use a target variance of 0
region      = 1;        % =0 for global (incl. Taiwan)
                        % =1 for California
                        % =2 for Japan
                        % =3 for China or Turkey
                        % =4 for Italy
Fault_Type  = 1;        % =0 for unspecified fault
                        % =1 for strike-slip fault
                        % =2 for normal fault
                        % =3 for reverse fault
                        
% Ground motion properties to require when selecting from the database. 
allowedVs30 = [-Inf Inf];     % upper and lower bound of allowable Vs30 values 
allowedMag  = [6 Inf];        % upper and lower bound of allowable magnitude values
allowedD    = [-Inf Inf];     % upper and lower bound of allowable distance values

% Miscellaneous other inputs
showPlots   = 1;        % =1 to plot results, =0 to suppress plots
seedValue   = 1;        % =0 for random seed in when simulating 
                        % response spectra for initial matching, 
                        % otherwise the specifed seedValue is used.
nTrials     = 20;       % number of iterations of the initial spectral 
                        % simulation step to perform
outputFile  = 'Output_File.dat'; % File name of the output file

% User inputs end here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the user-chosen database and format according to type of selection
% load the specified database
load(databaseFile) 

% Format appropriate ground motion metadata variables for single or two-
% component selection. Additional metadata is available in the databases
% and can be added here if desired 
% Note: These lines should be modified if using the BBP_EXSIM_meta_data.mat
% database file. See documentation for more details.
if arb == 1
    Filename    = [Filename_1; Filename_2];
    SaKnown     = [Sa_1; Sa_2]; 
    soil_Vs30   = [soil_Vs30; soil_Vs30]; 
    magnitude   = [magnitude; magnitude]; 
    closest_D   = [closest_D; closest_D]; 
    dirLocation = [dirLocation; dirLocation];
else % two-component selection
    Filename    = Filename_1;
    if RotD == 50 && exist('Sa_RotD50')
        SaKnown     = Sa_RotD50;
    elseif RotD == 100 && exist('Sa_RotD100')
        SaKnown     = Sa_RotD100;
    else
        fprintf(['Error--RotD' num2str(RotD) ' not provided in database \n\n'])
        % if data corresponding to user input RotD value does not exist,
        % use the geometric mean of two single-component directions
        SaKnown = sqrt(Sa_1.*Sa_2);
    end
end

%% Arrange available spectra in usable format and check for invalid values
% Create variable for known periods
knownPer = Periods; 

% Modify TgtPer to include T1 if running a conditional selection
if optInputs.cond == 1 && ~any(optInputs.TgtPer == optInputs.T1)
    optInputs.TgtPer = sort([optInputs.TgtPer optInputs.T1]);
end

% Match periods (known periods and target periods for error computations) 
% save the indicies of the matched periods in knownPer
indPer = zeros(length(optInputs.TgtPer),1);
for i=1:length(optInputs.TgtPer)
    [~ , indPer(i)] = min(abs(knownPer - optInputs.TgtPer(i)));
end

% Remove any repeated values from TgtPer and redefine TgtPer as periods 
% provided in databases
indPer = unique(indPer);
optInputs.TgtPer = knownPer(indPer);

% Identify the index of T1 within TgtPer 
[~, optInputs.indT1] = min(abs(optInputs.TgtPer - optInputs.T1));

%% Screen the records to be considered
recValidSa = ~all(SaKnown == -999,2); % remove invalid inputs
recValidSoil = soil_Vs30 > allowedVs30(1) & soil_Vs30 < allowedVs30(2);
recValidMag =  magnitude > allowedMag(1)  & magnitude < allowedMag(2);
recValidDist = closest_D > allowedD(1)    & closest_D < allowedD(2);

% find indicies of allowable records that will be searched
allowedIndex = find(recValidSoil & recValidMag & recValidDist & recValidSa); 

% Process available spectra
SaKnown = SaKnown(allowedIndex,:);       % allowed spectra defined at all periods
IMs.sampleBig = log(SaKnown(:,indPer));  % logarithmic spectral accelerations at target periods
optInputs.nBig = size(IMs.sampleBig,1);  % number of allowed spectra

fprintf('Number of allowed ground motions = %i \n \n', optInputs.nBig)
assert(optInputs.nBig >= optInputs.nGM, 'Warning: there are not enough allowable ground motions');

%% Compute target means and covariances of spectral values 
% Input the following variables and comment out ComputeTargets() for a
% user-defined target response spectrum 
% knownMeanReq = []; % (log) target mean response spectrum with length of knownPer
% knownCovReq = [];  % target covariance matrix with size [length(knownPer) length(knownPer)]

% Compute the mean response spectrum and target covariance matrix at all
% available periods of the database
[knownMeanReq, knownCovReq] = ComputeTargets(RotD, arb, indPer, knownPer, useVar, eps_bar, optInputs, ...
                                                M_bar, Rjb, Fault_Type, region, z1, Vs30);
                                            
% Define target spectrum and mean covariance matrix at target periods
% defined above
Tgts.meanReq = knownMeanReq(indPer); 
Tgts.covReq  = knownCovReq(indPer,indPer);

%% Simulate response spectra matching the above targets, for use in record selection
% Set initial seed for simulation
if seedValue ~= 0
    rng(seedValue); 
else
    rng('shuffle');
end

% Generate simulated response spectra with best matches to the target values
devTotalSim = zeros(nTrials,1);
for j=1:nTrials
    spectraSample{j} = exp(lhsnorm(Tgts.meanReq,Tgts.covReq,optInputs.nGM));
    sampleMeanErr = mean(log(spectraSample{j})) - Tgts.meanReq; % how close is the mean of the spectra to the target
    sampleStdErr = std(log(spectraSample{j})) - sqrt(diag(Tgts.covReq))'; % how close is the standard dev. of the spectra to the target
    sampleSkewnessErr = skewness(log(spectraSample{j}),1); % how close is the skewness of the spectra to zero (i.e., the target)
    devTotalSim(j) = optInputs.weights(1) * sum(sampleMeanErr.^2) + ...
                     optInputs.weights(2) * sum(sampleStdErr.^2)+ ...
                     optInputs.weights(3) * sum(sampleSkewnessErr.^2); % combine the three error metrics
end
[~, bestSample] = min(devTotalSim); % find the simulated spectra that best match the targets 
simulatedSpectra = spectraSample{bestSample};

%% Find best matches to the simulated spectra from ground-motion database
% Define the spectral accleration at T1 that all ground motions will be scaled to
optInputs.lnSa1 = Tgts.meanReq(optInputs.indT1); 
if optInputs.cond == 1 
    scaleFacIndex = optInputs.indT1;
else
    scaleFacIndex = (1:length(optInputs.TgtPer))';
end

% Initialize vectors
optInputs.recID = zeros(optInputs.nGM,1);
IMs.sampleSmall = [];
finalScaleFac = ones(optInputs.nGM,1);

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
    err(optInputs.recID(1:i-1)) = 1000000; % exclude previously-selected ground motions 
    err(scaleFac > optInputs.maxScale) = 1000000; % exclude ground motions requiring too large of a scale factor
   
    % find minimum-error ground motion    
    [tmp, optInputs.recID(i)] = min(err);
    assert(tmp < 1000000, 'Warning: problem with simulated spectrum. No good matches found');
    finalScaleFac(i) = scaleFac(optInputs.recID(i)); % store scale factor
    IMs.sampleSmall = [IMs.sampleSmall; log(exp(IMs.sampleBig(optInputs.recID(i),:))*scaleFac(optInputs.recID(i)))]; % store scaled log spectrum
end

% Compute target means and standard deviations and the means 
% and standard deviations of the originally selected ground motions 
Tgts.stdevs = sqrt(diag(Tgts.covReq))';
origMeans = mean(log(SaKnown(optInputs.recID,:).*repmat(finalScaleFac,1,size(SaKnown,2))));
origStdevs= std(log(SaKnown(optInputs.recID,:).*repmat(finalScaleFac,1,size(SaKnown,2))));

% Compute maximum percent error of selection relative to target means and
% standard deviations (do not compute standard deviation error at T1 for
% conditional selection)
meanErr = max(abs(exp(origMeans(indPer))-exp(Tgts.meanReq))./exp(Tgts.meanReq))*100;
stdErr = max(abs(origStdevs(indPer ~= indPer(optInputs.indT1))-Tgts.stdevs(1:end ~= optInputs.indT1))./Tgts.stdevs(1:end ~= optInputs.indT1))*100;

% Display the original maximum error between the selected gm and the target
fprintf('End of simulation stage \n')
fprintf('Max (across periods) error in median = %3.1f percent \n', meanErr); 
fprintf('Max (across periods) error in standard deviation = %3.1f percent \n \n', stdErr); 

%% Further optimize the ground motion selection, if desired
% if selected motions do not yet meet tolerances, further optimize
if meanErr > optInputs.tol || stdErr > optInputs.tol 
    [IMs.sampleSmall, finalRecords, finalScaleFac] = GreedyOpt(optInputs, Tgts, IMs);  
    % [IMs.sampleSmall, finalRecords, finalScaleFactors] = GreedyOptPar(optInputs, Tgts, IMs); % a version of the optimization function that uses parallel processing
else % otherwise, skip greedy optimization
    display('Greedy optimization was skipped based on user input tolerance.');
    finalRecords = optInputs.recID;
end

%% Plot results, if desired
if (showPlots)
    plotResults;
end

%% Output results to a text file 
% Produce a tab-delimited file with selected ground motions and scale factors. 
% For instructions on downloading the time histories, see the documentation
% files for each database. 

fin = fopen(outputFile,'w');
% print header information
fprintf(fin, '%s \n \n', getTimeSeries{1}, getTimeSeries{2}, getTimeSeries{3});
if arb == 1
    fprintf(fin,'%s \t %s \t %s \t %s \t %s \n','Record Number','Record Sequence Number','Scale Factor','File Name','URL');
elseif arb == 2
    fprintf(fin,'%s \t %s \t %s \t %s \t %s \t %s \t %s \n','Record Number','Record Sequence Number','Scale Factor','File Name Dir. 1','File Name Dir. 2', 'URL 1', 'URL 2');
end

% print record data
for i = 1 : length(finalRecords)
    rec = allowedIndex(finalRecords(i));
    if arb == 1 
        fprintf(fin,'%d \t %d \t %6.2f \t %s \t %s \n',i,rec,finalScaleFac(i),char(Filename{rec}),[char(dirLocation{rec}) char(Filename{rec})]); % Print relevant outputs
    else 
        fprintf(fin,'%d \t %d \t %6.2f \t %s \t %s \t %s \t %s \n',i,rec,finalScaleFac(i),char(Filename_1{rec}),char(Filename_2{rec}),[char(dirLocation{rec}) char(Filename_1{rec})],[char(dirLocation{rec}) char(Filename_2{rec})]);
    end
end

fclose(fin);