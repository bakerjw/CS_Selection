% This code is used to select conditional and unconditional ground motions.
% Conditional selection chooses structure- and site- specific
% ground motions. The target means and covariances are obtained
% corresponding to a pre-defined target scenario earthquake, and are
% obtained based on the CMS method. For unconditional selection, the target 
% means and covariances are still corresponding to a pre-defined target 
% scenario earthquake, but are not conditioned at a specified period of 
% interest. The target means and variances can be
% modified by the user. The ground motion selection algorithm used here can 
% be based on a single arbitrary ground motion component or two-component
% ground motions based on a flag set by the user. The single arbitrary 
% ground motion component is used for two-dimensional structural models and
% two-component ground motions are used for three-dimensional structural 
% models, as described in Jayaram et al. (2011). 
%
% created by Nirmal Jayaram, Ting Lin, Jack W. Baker
% Last Updated: 7 June 2010 
%
% modified by Cynthia Lee and Jack Baker
% modified 4/23/2015
% modified 7/7/2015 (Spring & Summer 2015)
% modified 11/14/2015
% Last Updated: 
% 
% Referenced manuscripts:
%
% N. Jayaram, T. Lin and and Baker, J. W. (2011). A computationally
% efficient ground-motion selection algorithm for matching a target
% response spectrum mean and variance, Earthquake Spectra, 27(3), 797-815.
%
% Baker, J.W. (2011). The Conditional Mean Spectrum: A tool for ground
% motion selection, ASCE Journal of Structural Engineering, 137(3), 322-331.
%
% ENTER NEW REPORT TITLE

%% OUTPUT VARIABLES
%
% finalRecords      : Record numbers of selected records
% finalScaleFactors : Scale factors
%
% The final cell in this m file shows how to plot the selected spectra
% using this information.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variable definitions and initial set of user inputs
% Specify a database with a needed ground motion data. Provided databases
% include 'NGA_W1_meta_data.mat', 'NGA_W2_meta_data.mat' and
% 'GP_sim_meta_data.mat' Further documentation of these databases can be 
% found at 'WorkspaceDocumentation***.txt'.
%
% Variable definitions for loading data:
% databaseFile : filename of the target database. 
% cond         : 0 to run unconditional selection
%                1 to run conditional
% arb          : 1 for single-component selection and arbitrary component sigma
%                2 for two-component selection and average component sigma
% RotD         : 50 to use SaRotD50 data
%              : 100 to use SaRotD100 data
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Certain variables are stored as parts of data structures which will later
% be incorporated in the optimization portion of the code. Required user
% input values are indicated for the user. Some of these are labeled as
% advanced user inputs and are set to default values which the user may
% choose to change. Other variables described below are calculated within
% this script or other functions. The data structures are as follows:
%
% optInputs - input values needed to run the optimization function:
%            optType    : For greedy optimization, the user will input a 0
%                         to use the sum of squared errors approach to 
%                         optimize the selected spectra, or a 1 to use 
%                         D-statistic calculations from the KS-test
%            tol        : User input percent error tolerance to determine
%                         whether or not optimization can be skipped (this
%                         will only be used for SSE optimization)
%            nGM        : Number of ground motions to be selected 
%            PerTgt     : Periods at which the target spectrum needs to be
%                         computed (logarithmically spaced)
%            T1         : Period at which spectra should be scaled and 
%                         matched 
%            nBig       : The number of allowed spectra
%            recID      : This is a vector of index values for chosen
%                         spectra
%            rec        : This is the index of T1, the conditioning period
%            isScaled   : Should spectra be scaled before matching (1 -YES,
%                         0-NO). 
%            maxScale   : Maximum allowable scale factor. 
%            penalty    : If a penalty needs to be applied to avoid selecting
%                         spectra that have spectral acceleration values 
%                         beyond 3 sigma at any of the periods, set a value
%                         here. Use 0 otherwise.
%            weights    : [Weight for error in mean, Weight for error 
%                         in standard deviation] e.g., [1.0,1.0] - equal 
%                         weight for both errors (only needed for SSE
%                         optimization)
%            nLoop      : This is a meta-variable of the algorithm. The 
%                         greedy improvement technique does several passes 
%                         over the available data set and tries to improve 
%                         the selected records. This variable denotes the 
%                         number of passes. Recommended value: 2.
% 
% Tgts      - The target values (means and covariances) being matched
%            meanReq : Estimated target response spectrum means (vector of
%                      spectral values, one at each period)
%            covReq  : Matrix of response spectrum covariances
%            means   : The meanReq vector re-formatted to be equal to
%                      exp(meanReq) in order to plot the means 
%                      properly on a log scale
%            stdevs    : A vector of standard deviations at each period
% 
% IMs       - The intensity measure values (from SaKnown) chosen and the 
%             values available:
%            sampleSmall : The matrix of spectra chosen at any point; this
%                          matrix will be redefined after optimization so 
%                          that it will no longer be a part of the "IMs" 
%                          data structure
%            sampleBig   : The matrix of allowed spectra 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable definitions for more user inputs: 
% SaKnown   : (N*P matrix)
%             This is a matrix of Sa values at different periods (P) for
%             available ground-motion time histories (N).
%             Usually, the structure's fundamental period.
% Tmin      : Represents the shortest period at which the target spectrum 
%             needs to be computed
% Tmax      : The longest period at which the target spectrum needs to be
%             computed
% checkCorr : If 1, this runs a code that compares the correlation
%             structure of the selected ground motions to the correlations
%             published by Baker and Jayaram (2008).
% nTrials   : number of sets of response spectra that are simulated. The
%             best set (in terms of matching means, variances and skewness
%             is chosen as the seed). The user can also optionally rerun
%             this segment multiple times before deciding to proceed with
%             the rest of the algorithm. It is to be noted, however, that
%             the greedy improvement technique significantly improves the
%             match between the means and the variances subsequently.
% seedValue : For repeatability. For a particular seedValue not equal to
%             zero, the code will output the same set of ground motions.
%             The set will change when the seedValue changes. If set to
%             zero, the code randomizes the algorithm and different sets of
%             ground motions (satisfying the target mean and variance) are
%             generated each time.
% showPlots : 0 to suppress plots, 1 to show plots
% outputFile: File name of the output file
% perKnown  : The set of P periods.
% allowedVs30 : Only records with Vs30 values within this range will be
%               searched. For the simulated database, all Vs30 values are
%               863 m/s and this range must contain 863 m/s or the
%               algorithm will not run.
% allowedMag  : Only records with magnitudes within this range will be
%               searched
% allowedD    : Only records with closest distances within this range will
%               be searched
%
% If a database other than the NGA database is used, also define the
% following variables:
% 
% soil_Vs30        : Soil Vs30 values corresponding to all the records
% magnitude        : Magnitude of all the records
% closest_D        : Closest distance for all the records
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User inputs begin here
% Choose data set and type of selection the user should note that the
% original NGA database does not contain RotD100 values for two-component
% selection

databaseFile         = 'NGA_W2_meta_data';
optInputs.cond       = 1;
arb                  = 1; 
RotD                 = 50; 

% Choose number of ground motions and set requirements for periods
optInputs.nGM        = 20;
optInputs.T1         = 1.5; 
Tmin                 = 0.1;
Tmax                 = 10;

% Choose to show or suppress plots
showPlots            = 1;

% MORE user input in the next cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User inputs for determination of target mean and covariances
% The Campbell and Bozorgnia (2008) ground-motion model is used in this
% code. The input variables defined below are the inputs required for this
% model. The user can change the ground-motion model as long as any
% additional information that may be required by the new model is provided.
% Refer to Jayaram et al.(2011) and Baker (2011) for details on the 
% method used for obtaining the target means and covariances. 

% The code provides the user an option to not match the target variance.
% This is done by setting the target variance to zero so that each selected
% response spectrum matches the target mean response spectrum.

% Variable definitions

% M_bar     : Magnitude of the target scenario earthquake
% R_bar     : Distance corresponding to the target scenario earthquake
% eps_bar   : Epsilon value for the target scenario (for conditional
%             selection; for unconditional selection, leave as any value or
%             set equal to 0)
% Vs30      : Average shear wave velocity in the top 30m of the soil, used
%             to model local site effects (m/s)
% Ztor      : Depth to the top of coseismic rupture (km)
% delta     : Average dip of the rupture place (degree)
% lambda    : Rake angle (degree)
% Zvs       : Depth to the 2.5 km/s shear-wave velocity horizon (km)
% useVar    : 0 to set target variance to zero, 1 to compute target
%             variance using ground-motion model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs begin here

M_bar       = 7.4;
R_bar       = 10; 
eps_bar     = 2; % for conditional selection
Vs30        = 400;
Ztor        = 0;
delta       = 90;
lambda      = 180;
Zvs         = 2;
useVar      = 1;

Rrup        = R_bar; 
Rjb         = R_bar; 

%% Advanced user inputs  
% The definitions for these inputs are documented in the sections above.
% Most users will likely keep these default values.

% Choose limits to screen databases
allowedVs30          = [200 600]; 
allowedMag           = [5 7];
allowedD             = [0 30]; 

% Advanced user inputs for optimization 
optInputs.PerTgt     = logspace(log10(Tmin),log10(Tmax),30);
optInputs.isScaled   = 1;
optInputs.maxScale   = 4;
optInputs.tol        = 15; 
optInputs.weights    = [1.0 2.0];
optInputs.nLoop      = 2;
optInputs.penalty    = 0;
optInputs.optType    = 0; 

% Miscellaneous advanced inputs
seedValue            = 1; % default will be set to 0
nTrials              = 20;
checkCorr            = 1;
outputFile           = 'Output_File.dat';

% User inputs end here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the user-chosen database and format according to type of selection
% load the specified database
load(databaseFile) 

% Format appropriate variables for single or two-component selection
if arb == 1
    Filename    = [Filename_1; Filename_2];
    SaKnown     = [Sa_1; Sa_2]; 
    soil_Vs30   = [soil_Vs30; soil_Vs30]; 
    magnitude   = [magnitude; magnitude]; 
    closest_D   = [closest_D; closest_D]; 
else % two-component selection
    Filename    = Filename_1;
    if RotD == 50 && exist('Sa_RotD50')
        SaKnown     = Sa_RotD50;
    elseif RotD == 100 && exist('Sa_RotD100')
        SaKnown     = Sa_RotD100;
    else
        fprintf(['Error--RotD' num2str(RotD) ' not provided in database \n\n'])
    end
end

% Create variable for known periods
perKnown = Periods;

% More fields available in databases that can also be used in screening 
% (e.g. the ones shown below)

% mechanism = [mechanism; mechanism];
% lowest_usable_freq = [lowest_usable_freq; lowest_usable_freq];
% distance_jb = [distance_jb; distance_jb];
% distance_hyp = [distance_hyp; distance_hyp];
% distance_epi = [distance_epi; distance_epi];
% distance_campbell = [distance_campbell; distance_campbell];

%% Arrange available spectra in usable format and check for invalid values

% Modify PerTgt to include T1 if running a conditional selection
if optInputs.cond == 1 && ~any(optInputs.PerTgt == optInputs.T1)
    optInputs.PerTgt = [optInputs.PerTgt(optInputs.PerTgt<optInputs.T1)...
        optInputs.T1 optInputs.PerTgt(optInputs.PerTgt>optInputs.T1)];
end

% Match periods (known periods and periods for error computations) save the
% indicies of the matched periods in perKnown
recPer = zeros(length(optInputs.PerTgt),1);
for i=1:length(optInputs.PerTgt)
    [~ , recPer(i)] = min(abs(perKnown - optInputs.PerTgt(i)));
end

% Remove any repeated values from PerTgt and redefine PerTgt as periods 
% provided in databases
recPer = unique(recPer);
optInputs.PerTgt = perKnown(recPer);

% Identify the index of T1 within PerTgt and the final number of periods in
% PerTgt
[~, optInputs.rec] = min(abs(optInputs.PerTgt - optInputs.T1));
numPer = length(optInputs.PerTgt);

% Screen the records to be considered
recValidSa = ~all(SaKnown == -999,2); % remove invalid inputs
recValidSoil = soil_Vs30 > allowedVs30(1) & soil_Vs30 < allowedVs30(2);
recValidMag = magnitude > allowedMag(1) & magnitude < allowedMag(2);
recValidDist = closest_D > allowedD(1) & closest_D < allowedD(2);

% only the allowable records will be searched
allowedIndex = find(recValidSoil & recValidMag & recValidDist & recValidSa); 

SaKnown = SaKnown(allowedIndex,:);
IMs.sampleBig = SaKnown(:,recPer);

% Processing available spectra
IMs.sampleBig = log(IMs.sampleBig);
optInputs.nBig = size(IMs.sampleBig,1);

fprintf('Number of allowed ground motions = %i \n \n', length(allowedIndex))


%% Determine target spectra and database correlations using ground-motion model 

% arrange periods for which correlations will be calculated
perKnownCorr = perKnown;
if optInputs.cond == 1 && ~any(perKnown == optInputs.T1)
    perKnownCorr = [perKnown(perKnown<optInputs.T1) optInputs.T1 perKnown(perKnown>optInputs.T1)];
end
perKnownCorr = perKnownCorr(perKnownCorr <=10);

% compute the response spectrum values from the Campbell and
% Bozorgnia GMPE at all periods for which correlations will be calculated
[saCorr, sigmaCorr] = CB_2008_nga (M_bar, perKnownCorr, Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb);

% compute the target means, covariances, and correlations 
% return target means and covariances within the structure "Tgts" 
[scaleFacIndex, corrReq, Tgts, optInputs] = ComputeTargets(recPer, perKnown, perKnownCorr, saCorr,...
                                                sigmaCorr, useVar, eps_bar, optInputs);
                                            
%% Simulate response spectra using Monte Carlo Random Simulation/Latin Hypercube Sampling

% Set initial seed for simulation
if seedValue ~= 0
    rng(seedValue); 
else
    rng('shuffle');
end

% Generate simulated response spectra with best matches to the target values
devTotalSim = zeros(nTrials,1);
for j=1:nTrials
    gmCell{j} = zeros(optInputs.nGM,length(optInputs.PerTgt));
    gmCell{j} = exp(lhsnorm(Tgts.meanReq,Tgts.covReq,optInputs.nGM));
    devMeanSim = mean(log(gmCell{j})) - Tgts.meanReq;
    devSkewSim = skewness(log(gmCell{j}),1);
    devSigSim = std(log(gmCell{j})) - sqrt(diag(Tgts.covReq))';
    devTotalSim(j) = optInputs.weights(1) * sum(devMeanSim.^2) + ...
                     optInputs.weights(2) * sum(devSigSim.^2)+ 0.1 * ...
                     (optInputs.weights(1)+optInputs.weights(2)) * sum(devSkewSim.^2);
end
[tmp, recUse] = min(abs(devTotalSim));
gm = gmCell{recUse};

%% Find best matches to the simulated spectra from ground-motion database

% Initialize vectors
optInputs.recID = zeros(optInputs.nGM,1);
IMs.sampleSmall = [];
finalScaleFac = ones(optInputs.nGM,1);

% Match spectra from database to simulated spectra by calculating
% difference between them. If a scale factor is greater than maximum scale
% factor or if spectra is already chosen, set error to 1000000
for i = 1:optInputs.nGM
    err = zeros(optInputs.nBig,1);
    scaleFac = ones(optInputs.nBig,1);
    
    for j=1:optInputs.nBig
        if (any(optInputs.recID == j)) 
            err(j) = 1000000;
        elseif optInputs.isScaled == 1
            scaleFac(j) = sum(exp(IMs.sampleBig(j,scaleFacIndex)).*gm(i,scaleFacIndex))/sum(exp(IMs.sampleBig(j,scaleFacIndex)).^2);
            if scaleFac(j) > optInputs.maxScale
                err(j) = 1000000;
            else
                err(j) = sum((log(exp(IMs.sampleBig(j,:))*scaleFac(j)) - log(gm(i,:))).^2); 
            end
        else
            err(j) = sum((IMs.sampleBig(j,:) - log(gm(i,:))).^2);
        end
    end
    
    [tmp, optInputs.recID(i)] = min(err);
    if tmp >= 1000000
       display('Warning: Possible problem with simulated spectrum. No good matches found');
       display(optInputs.recID(i));
    end
    if (optInputs.isScaled == 1)
       finalScaleFac(i) = scaleFac(optInputs.recID(i));
    else
       finalScaleFac(i) = 1;
    end
    IMs.sampleSmall = [IMs.sampleSmall;log(exp(IMs.sampleBig(optInputs.recID(i),:))*scaleFac(optInputs.recID(i)))];

end

%% Skip greedy optimization if the user-defined tolerance for maximum 
% percent error between the target and sample means and standard deviations
% have been attained

% Format and store the target means and standard deviations and the means 
% and standard deviations of the originally selected ground motions 
Tgts.means = exp(Tgts.meanReq);
Tgts.stdevs = sqrt(diag(Tgts.covReq))';
origMeans = exp(mean(IMs.sampleSmall));
origStdevs = std(IMs.sampleSmall);

% Define all other periods besides T1
notT1 = find(optInputs.PerTgt ~= optInputs.PerTgt(optInputs.rec)); 

% Compute maximum percent error from target
meanErr = max(abs(origMeans-Tgts.means)./Tgts.means)*100;
stdErr = max(abs(origStdevs(notT1)-Tgts.stdevs(notT1))./Tgts.stdevs(notT1))*100;

% Display the original maximum error between the selected gm and the target
fprintf('End of simulation stage \n')
fprintf('Max (across periods) error in median = %3.1f percent \n', meanErr); 
fprintf('Max (across periods) error in standard deviation = %3.1f percent \n \n', stdErr); 

%% Greedy subset modification procedure
% Call the optimization function if the user defined tolerance has not been
% reached. Each nLoop, the error is recalculated and the loop will break if
% the error has been reached

% numWorkers = 2;
% parobj = parpool(numWorkers);
if meanErr > optInputs.tol || stdErr > optInputs.tol 
    [sampleSmall, finalRecords, finalScaleFactors] = GreedyOpt(optInputs, Tgts, IMs);
    IMs.sampleSmall = sampleSmall;
    
else % otherwise, skip greedy optimization
    display('Greedy optimization was skipped based on user input tolerance.');
    finalRecords = optInputs.recID;
    finalScaleFactors = finalScaleFac;
end

% delete(parobj);

%% Spectra Plots

if (showPlots)
    plotResults;
end

%% Output data to file (best viewed with a text editor)
% To obtain the time histories for the ground motions represented in each
% database, proceed to the links provided in the output file or follow the
% instructions provided at the top of the file. For more specific
% instructions for downloading the time histories, see the documentation
% files for each database. 

fin = fopen(outputFile,'w');
% print header information
fprintf(fin, '%s \n \n', getTimeSeries{1}, getTimeSeries{2}, getTimeSeries{3});
if arb == 1
    fprintf(fin,'%s \t %s \t %s \t %s \n','Record Number','Scale Factor','File Name','URL');
elseif arb == 2
    fprintf(fin,'%s \t %s \t %s \t %s \t %s \t %s \t %s \n','Record Number','NGA Record Sequence Number','Scale Factor','File Name Dir. 1','File Name Dir. 2','URL Dir 1','URL Dir 2');
end

% print record data
for i = 1 : length(finalRecords)
    rec = allowedIndex(finalRecords(i));
    if arb == 1
        fprintf(fin,'%d \t %6.2f \t %s \t %s \n',i,finalScaleFactors(i),Filename{rec},[dirLocation Filename{rec}]); % Print relevant outputs
    else 
        fprintf(fin,'%d \t %d \t %6.2f \t %s \t %s \t %s \t %s \n',i,rec,finalScaleFactors(i),Filename_1{rec},Filename_2{rec},[dirLocation Filename_1{rec}],[dirLocation Filename_2{rec}]);
    end
end

fclose(fin);