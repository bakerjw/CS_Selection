function [ SaKnown, selectionParams, indPer, knownPer, Filename, dirLocation, getTimeSeries, allowedIndex ] = screen_database(selectionParams, databaseFile, allowedRecs )
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
if selectionParams.arb == 1 % single-component selection -- treat each component as a seaparate candidate
    Filename    = [Filename_1; Filename_2];
    SaKnown     = [Sa_1; Sa_2]; 
    soil_Vs30   = [soil_Vs30; soil_Vs30]; 
    magnitude   = [magnitude; magnitude]; 
    closest_D   = [closest_D; closest_D]; 
    dirLocation = [dirLocation; dirLocation];
else % two-component selection
    Filename    = [Filename_1 Filename_2];
    if selectionParams.RotD == 50 && exist('Sa_RotD50')
        SaKnown     = Sa_RotD50;
    elseif selectionParams.RotD == 100 && exist('Sa_RotD100')
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
if selectionParams.cond == 1 && ~any(selectionParams.TgtPer == selectionParams.T1)
    selectionParams.TgtPer = sort([selectionParams.TgtPer selectionParams.T1]);
end

% Match periods (known periods and target periods for error computations) 
% save the indicies of the matched periods in knownPer
indPer = zeros(length(selectionParams.TgtPer),1);
for i=1:length(selectionParams.TgtPer)
    [~ , indPer(i)] = min(abs(knownPer - selectionParams.TgtPer(i)));
end

% Remove any repeated values from TgtPer and redefine TgtPer as periods 
% provided in databases
indPer = unique(indPer);
selectionParams.TgtPer = knownPer(indPer);

% Identify the index of T1 within TgtPer 
[~, selectionParams.indT1] = min(abs(selectionParams.TgtPer - selectionParams.T1));

%% Screen the records to be considered
recValidSa = ~all(SaKnown == -999,2); % remove invalid inputs
recValidSoil = soil_Vs30 > allowedRecs.Vs30(1) & soil_Vs30 < allowedRecs.Vs30(2);
recValidMag =  magnitude > allowedRecs.Mag(1)  & magnitude < allowedRecs.Mag(2);
recValidDist = closest_D > allowedRecs.D(1)    & closest_D < allowedRecs.D(2);

% flag indicies of allowable records that will be searched
allowedIndex = find(recValidSoil & recValidMag & recValidDist & recValidSa); 

% resize SaKnown to include only allowed records
SaKnown = SaKnown(allowedIndex,:);       

% count number of allowed spectra
selectionParams.nBig = length(allowedIndex);  

fprintf('Number of allowed ground motions = %i \n \n', selectionParams.nBig)
assert(selectionParams.nBig >= selectionParams.nGM, 'Warning: there are not enough allowable ground motions');


end

