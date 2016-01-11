function [ sampleSmall, finalRecords, finalScaleFac ] = GreedyOpt( optInputs, Tgts, IMs )
% This function will perform a greedy optimization on a set of ground
% motions using the sum of squared errors approach to check the set of
% selected ground motions against target means and variances
%   The inputs for the function are three structures that are also defined 
%   in the Select_Ground_Motions script:
%       optInputs - 
%                  cond       : 0 to run unconditional selection; 1 to run 
%                               conditional
%                  nGM        : Number of ground motions to be selected
%                  T1         : Period at which spectra should be scaled and matched
%                  TgtPer     : Periods at which the target spectrum needs to 
%                               be computed (logarithmically spaced)
%                  indT1      : Index of T1 within PerTgt
%                  lnSa1      : Spectral accleration at T1 that all ground 
%                               motions will be scaled to
%                  nBig       : The number of spectra that will be searched
%                  isScaled   : The user will input 1 to allow records to be 
%                               scaled, and input 0 otherwise 
%                  maxScale   : The maximum allowable scale factor
%                  indT1      : This is the index of T1, the conditioning period
%                  recID      : This is a vector of index values for chosen
%                               spectra
%                  tol        : User input percent error tolerance to determine
%                               whether or not optimization can be skipped (this
%                               will only be used for SSE optimization)
%                  optType    : For greedy optimization, the user will input a 0
%                               to use the sum of squared errors approach to 
%                               optimize the selected spectra, or a 1 to use 
%                               D-statistic calculations from the KS-test
%                  penalty    : If a penalty needs to be applied to avoid selecting
%                               spectra that have spectral acceleration values 
%                               beyond 3 sigma at any of the periods, set a value
%                               here. Use 0 otherwise.
%                  weights    : [Weight for error in mean, Weight for error 
%                              in standard deviation] e.g., [1.0,1.0] (used only
%                              in SSE optimization)
%                  nLoop      : Number of loops of optimizations to perform.
%                              Default value = 2
%
%       Tgts      - The target values (means and covariances) being matched
%                  meanReq    : Estimated target response spectrum means (vector of
%                               spectral values, one at each period)
%                  covReq     : Matrix of response spectrum covariances
%                  means      : The meanReq vector re-formatted to be equal to
%                               exp(meanReq) in order to plot the means 
%                               properly on a log scale
%                  stdevs     : A vector of standard deviations at each period
% 
%       IMs       - The intensity measure values (from SaKnown) chosen and the 
%                   values available:
%                  sampleSmall : matrix of selected response spectra 
%                  sampleBig   : The matrix of spectra that will be searched
% 

% Redefine sampleSmall since the size will be continuously changing as the
% optimization runs. The output sampleSmall will not be returned as part of
% a data structure, but as a new and separate variable
sampleSmall = IMs.sampleSmall;

display('Please wait...This algorithm takes a few minutes depending on the number of records to be selected');
if optInputs.cond == 0
    display('The algorithm is slower when scaling is used');
end
if optInputs.optType == 1
    display('The algorithm is slower when optimizing with the KS-test Dn statistic');
end

% if optimizing the ground motions by calculating the Dn value, first
% calculate the emperical CDF values (which will be the same at each
% period) 
if optInputs.optType == 1
    emp_cdf = linspace(0,1,optInputs.nGM+1);
end

% Initialize scale factor vector for each iteration and a scale factor
% vector for the final selected records 
scaleFac = ones(optInputs.nBig,1);
finalScaleFac = ones(optInputs.nGM,1);

for k=1:optInputs.nLoop % Number of passes
    for i=1:optInputs.nGM % Selects nGM ground motions
        display([num2str(round(((k-1)*optInputs.nGM + i-1)/(optInputs.nLoop*optInputs.nGM)*100)) '% done']);
        devTotal = zeros(optInputs.nBig,1); % initialize measure of deviation for each selected ground motion
        sampleSmall(i,:) = []; % remove initially selected record to be replaced
        optInputs.recID(i,:) = []; 
        
        if optInputs.isScaled == 1
            if optInputs.cond == 1
                scaleFac = exp(optInputs.lnSa1)./exp(IMs.sampleBig(:,optInputs.indT1));
            elseif optInputs.cond == 0
                [scaleFac, devTotal] = bestScaleFactor(IMs.sampleBig, sampleSmall, Tgts.meanReq, Tgts.stdevs, optInputs.weights, optInputs.maxScale);
            end
        end
        % Try to add a new spectra to the subset list
        for j = 1:optInputs.nBig
            sampleSmall = [sampleSmall;IMs.sampleBig(j,:)+log(scaleFac(j))];
            
            % Calculate the appropriate measure of deviation and store in
            % devTotal (the D-statistic or the combination of mean and
            % sigma deviations)
            if optInputs.optType == 0
                if optInputs.cond == 1 | (optInputs.cond == 0 && optInputs.isScaled == 0)
                    % Compute deviations from target
                    devMean = mean(sampleSmall) - Tgts.meanReq;
                    devSig = std(sampleSmall) - Tgts.stdevs;
                    devTotal(j) = optInputs.weights(1) * sum(devMean.^2) + optInputs.weights(2) * sum(devSig.^2);
                end
                % Penalize bad spectra (set penalty to zero if this is not required)
                if optInputs.penalty ~= 0
                    for m=1:size(sampleSmall,1)
                        devTotal(j) = devTotal(j) + sum(abs(exp(sampleSmall(m,:))>exp(Tgts.meanReq+3*Tgts.stdevs'))) * optInputs.penalty;
                    end
                end
                
            elseif optInputs.optType == 1   
                sortedlnSa = [min(sampleSmall); sort(sampleSmall)];
                norm_cdf = normcdf(sortedlnSa,repmat(Tgts.meanReq,optInputs.nGM+1,1),repmat(Tgts.stdevs,optInputs.nGM+1,1));
                Dn = max(abs(repmat(emp_cdf',1,length(optInputs.PerTgt)) - norm_cdf));
                devTotal(j) = sum(Dn);   
            end
            
            % Scale factors for either type of optimization should not
            % exceed the maximum
            if (scaleFac(j) > optInputs.maxScale)
                devTotal(j) = 1000000;
            end
            
            % Should cause improvement and record should not be repeated
            if (any(optInputs.recID == j))
                devTotal(j) = 1000000;
            end
            sampleSmall = sampleSmall(1:end-1,:);
        end
        
        [~ , minID] = min(devTotal);
        % Add new element in the right slot
        if optInputs.isScaled == 1
            finalScaleFac(i) = scaleFac(minID);
        else
            finalScaleFac(i) = 1;
        end
        sampleSmall = [sampleSmall(1:i-1,:);IMs.sampleBig(minID,:)+log(scaleFac(minID));sampleSmall(i:end,:)];
        optInputs.recID = [optInputs.recID(1:i-1);minID;optInputs.recID(i:end)];
        
    end
    
    % Can the optimization be stopped after this loop based on the user
    % specified tolerance? Recalculate new standard deviations of new
    % sampleSmall and then recalculate new maximum percent errors of means
    % and standard deviations 
    if optInputs.optType == 0
        stdevs = std(sampleSmall);
        meanErr = max(abs(exp(mean(sampleSmall))-Tgts.means)./Tgts.means)*100;
        stdErr = max(abs(stdevs(1:end ~= optInputs.indT1) - Tgts.stdevs(1:end ~= optInputs.indT1))./Tgts.stdevs(1:end ~= optInputs.indT1))*100;
        fprintf('Max (across periods) error in median = %3.1f percent \n', meanErr); 
        fprintf('Max (across periods) error in standard deviation = %3.1f percent \n \n', stdErr);
        
        % If error is now within the tolerance, break out of the
        % optimization
        if meanErr < optInputs.tol && stdErr < optInputs.tol
            display('The percent errors between chosen and target spectra are now within the required tolerances.');
            break;
        end
    end
    
    fprintf('End of loop %1.0f of %1.0f \n', k, optInputs.nLoop) 
end

display('100% done');

% Output information
finalRecords = optInputs.recID;
finalScaleFac = finalScaleFac';




end

