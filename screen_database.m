function [ SaKnown, selectionParams, indPer, knownPer, metadata ] = screen_database(selectionParams, allowedRecs )
% Load a database of ground motion data, and screen it to identify usable
% ground motions for potential selection


% load the specified database
load(['Databases/' selectionParams.databaseFile]) 
metadata.getTimeSeries = getTimeSeries;
metadata.dirLocation = dirLocation;

% Format appropriate ground motion metadata variables for single or two-
% component selection. Additional metadata is available in the databases
% and can be added here if desired 
% Note: These lines should be modified if using the BBP_EXSIM_meta_data.mat
% database file. See documentation for more details.
if selectionParams.arb == 1 && selectionParams.matchV ~= 1 % single-component selection -- treat each component as a seaparate candidate
    if isempty(Sa_2) % if 2nd component doesn't exist (e.g., EXSIM data)
        metadata.Filename    = Filename_1;
        metadata.compNum     = [ones(size(magnitude))];
        metadata.dirLocation = dirLocation;
        metadata.recNum      = [1:length(magnitude)]';
        SaKnown              = Sa_1;
    else % 2nd component exists
        metadata.Filename    = [Filename_1; Filename_2];
        metadata.compNum     = [ones(size(magnitude)); 2*ones(size(magnitude))];
        metadata.recNum      = [1:length(magnitude) 1:length(magnitude)]';
        SaKnown     = [Sa_1; Sa_2];
        soil_Vs30   = [soil_Vs30; soil_Vs30];
        magnitude   = [magnitude; magnitude];
        closest_D   = [closest_D; closest_D];
        NGA_num     = [NGA_num; NGA_num];
        metadata.dirLocation = [dirLocation; dirLocation];
    end
elseif selectionParams.arb == 2 && selectionParams.matchV ~= 1 % two-component selection
    metadata.Filename    = [Filename_1 Filename_2];
    metadata.dirLocation = dirLocation;
    metadata.recNum      = [1:length(magnitude)]';
    if selectionParams.RotD == 50 && exist('Sa_RotD50')
        SaKnown     = Sa_RotD50;
    elseif selectionParams.RotD == 100 && exist('Sa_RotD100')
        SaKnown     = Sa_RotD100;
    else
        display(['Error--RotD' num2str(RotD) ' not provided in database'])
        % If data corresponding to user input RotD value does not exist,
        % optionally use the geometric mean of the single-component Sa's:
        % SaKnown = sqrt(Sa_1.*Sa_2);
    end
else
    metadata.Filename    = [Filename_1 Filename_2 Filename_vert];
    metadata.dirLocation = dirLocation;
    metadata.recNum      = [1:length(magnitude)];
    if selectionParams.RotD == 50 && exist('Sa_RotD50')
        SaKnown     = Sa_RotD50;
    elseif selectionParams.RotD == 100 && exist('Sa_RotD100')
        SaKnown     = Sa_RotD100;
    else
        display(['Error--RotD' num2str(RotD) ' not provided in database'])
        % If data corresponding to user input RotD value does not exist,
        % optionally use the geometric mean of the single-component Sa's:
        % SaKnown = sqrt(Sa_1.*Sa_2);
    end
    SaKnownV     = Sa_vert;
end

%% Arrange available spectra in usable format and check for invalid values

% Create variable for known periods
idxPer = find(Periods <= 10); % throw out periods > 10s, as they cause problems with the GMPE evaluation
knownPer = Periods(idxPer); 

% Modify TgtPer to include Tcond if running a conditional selection
if selectionParams.cond == 1 && ~any(selectionParams.TgtPer == selectionParams.Tcond)
    selectionParams.TgtPer = sort([selectionParams.TgtPer selectionParams.Tcond]);
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

% Identify the index of Tcond within TgtPer 
[~, selectionParams.indTcond] = min(abs(selectionParams.TgtPer - selectionParams.Tcond));

% If selecting V components, revise TgtPerV
if selectionParams.matchV == 1 
    if selectionParams.cond == 1 && ~any(selectionParams.TgtPerV == selectionParams.Tcond)
        selectionParams.TgtPerV = sort([selectionParams.TgtPerV selectionParams.Tcond]);
    end
    indPerV = zeros(length(selectionParams.TgtPerV),1);
    for i=1:length(selectionParams.TgtPerV)
        [~ , indPerV(i)] = min(abs(knownPer - selectionParams.TgtPerV(i)));
    end
    indPerV = unique(indPerV);
    selectionParams.TgtPerV = knownPer(indPerV);
end

%% Screen the records to be considered
% If selecting V components, ensure that each of remaining records contains all 3 components
if selectionParams.matchV == 1
    recValidSa = ~any(Sa_1==-999,2) & ~any(Sa_2==-999,2) & ~any(ismember(Sa_vert,[-999 inf]),2);
else
    recValidSa = ~all(SaKnown == -999,2); % remove invalid inputs
end
recValidSoil = soil_Vs30 > allowedRecs.Vs30(1) & soil_Vs30 < allowedRecs.Vs30(2);
recValidMag =  magnitude > allowedRecs.Mag(1)  & magnitude < allowedRecs.Mag(2);
recValidDist = closest_D > allowedRecs.D(1)    & closest_D < allowedRecs.D(2);

recValidIdx = ~ismember(metadata.recNum,allowedRecs.idxInvalid);

% flag indicies of allowable records that will be searched
metadata.allowedIndex = find(recValidSoil & recValidMag & recValidDist & recValidSa & recValidIdx); 

% resize SaKnown to include only allowed records
SaKnown = SaKnown(metadata.allowedIndex,idxPer);

% If selecting V components, save new variables in selectionParams
if selectionParams.matchV == 1
    SaKnownV = SaKnownV(metadata.allowedIndex,idxPer);
    selectionParams.SaKnownV = SaKnownV;
    selectionParams.indPerV = indPerV;
end

% count number of allowed spectra
selectionParams.nBig = length(metadata.allowedIndex);  

display(['Number of allowed ground motions = ' num2str(selectionParams.nBig)])
assert(selectionParams.nBig >= selectionParams.nGM, 'Warning: there are not enough allowable ground motions');


end

