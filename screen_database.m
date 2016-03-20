function [ SaKnown, optInputs, indPer, knownPer, Filename, dirLocation, getTimeSeries, allowedIndex ] = screen_database(optInputs, databaseFile, arb, RotD, allowedVs30, allowedMag, allowedD )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Load the ground motion database and set up data matrices

% load the specified database
load(['Databases/' databaseFile]) 

% Format appropriate ground motion metadata variables for single or two-
% component selection. Additional metadata is available in the databases
% and can be added here if desired 
% Note: These lines should be modified if using the BBP_EXSIM_meta_data.mat
% database file. See documentation for more details.
if arb == 1 % single-component selection -- treat each component as a seaparate candidate
    Filename    = [Filename_1; Filename_2];
    SaKnown     = [Sa_1; Sa_2]; 
    soil_Vs30   = [soil_Vs30; soil_Vs30]; 
    magnitude   = [magnitude; magnitude]; 
    closest_D   = [closest_D; closest_D]; 
    dirLocation = [dirLocation; dirLocation];
else % two-component selection
    Filename    = [Filename_1 Filename_2];
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

% flag indicies of allowable records that will be searched
allowedIndex = find(recValidSoil & recValidMag & recValidDist & recValidSa); 


end

