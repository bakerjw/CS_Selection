function [ sampleSmall, finalRecords, finalScaleFactors ] = GreedyOpt( optInputs, Tgts, IMs )
% This function will perform a greedy optimization on a set of ground
% motions using the sum of squared errors approach to check the set of
% selected ground motions against target means and variances
%   The inputs for the function are three structures that are also defined 
%   in the Select_Ground_Motions script:
%       optInputs   : There are twelve general inputs needed for
%                     optimization -
%                           nLoop    : number of times full optimization is
%                                      run                       
%                           nGM      : number of ground motions to select
%                           cond     : 0 for conditional selection, else 1
%                           T1       : conditonal selection - period to
%                                      match
%                           isScaled : to scale spectra or not
%                           penalty  : add a penalty for bad spectra
%                           weights  : weight for error in mean and std dev
%                           PerTgt   : periods at which target spectra is
%                                      calculated
%                           maxScale : maximum scaling factor for spectra
%                           nBig     : number of allowed spectra
%                           recID    : indicies of 
%                           tol      : user-defined percent error tolerance
%                           optType  : 0 for sum of squared errors
%                                      calculations, 1 for KS-test 
%                                      Dn value calculations
%
%       Tgts        : This structure represents the target values the user
%                     would like to match -
%                           meanReq     : The target means (vector)
%                           covReq      : The target covariances (matrix)
%                           means       : exp(meanReq) - formatted for
%                                         plotting
%                           sigs        : sqrt(diag(covReq)) - formatted
%                                         for plotting
%       IMs         : This structure represents the intensity measures
%                     being optimized -
%                           sampleSmall   : The matrix of ground motions
%                                           being optimized
%                           sampleBig     : The matrix of available ground
%                                           motions the set is being chosen 
%                                           from

% Redefine sampleSmall since the size will be continuously changing as the
% optimization runs. The output sampleSmall will not be returned as part of
% a data structure, but as a new and separate variable
sampleSmall = IMs.sampleSmall;

% if optimizing the ground motions by calculating the Dn value, first
% calculate the emperical CDF values and initialize a vector of Dn values
if optInputs.optType == 1
    emp_cdf = linspace(0,1,optInputs.nGM+1);
    Dn = zeros(length(optInputs.PerTgt),1);
end

display('Please wait...This algorithm takes a few minutes depending on the number of records to be selected');
if optInputs.cond == 0
    display('The algorithm is slower when scaling is used');
end

for k=1:optInputs.nLoop % Number of passes
    
    for i=1:optInputs.nGM % Selects nGM ground motions
        display([num2str(round(((k-1)*optInputs.nGM + i-1)/(optInputs.nLoop*optInputs.nGM)*100)) '% done']);
        
        devTotal = zeros(optInputs.nBig,1);
        sampleSmall(i,:) = [];
        optInputs.recID(i,:) = [];
        
        if optInputs.isScaled == 1
            if optInputs.cond == 1
                scaleFac = exp(optInputs.lnSa1)./exp(IMs.sampleBig(:,optInputs.PerTgt == optInputs.T1));
            elseif optInputs.cond == 0
                [scaleFac, devTotal] = bestScaleFactor(IMs.sampleBig, sampleSmall, Tgts.meanReq, Tgts.sigs, optInputs.weights, optInputs.maxScale);
            end
        else
            scaleFac = ones(optInputs.nBig,1);
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
                    devSig = std(sampleSmall) - Tgts.sigs;
                    devTotal(j) = optInputs.weights(1) * sum(devMean.^2) + optInputs.weights(2) * sum(devSig.^2);
                end
                % Penalize bad spectra (set penalty to zero if this is not required)
                if optInputs.penalty ~= 0
                    for m=1:size(sampleSmall,1)
                        devTotal(j) = devTotal(j) + sum(abs(exp(sampleSmall(m,:))>exp(Tgts.meanReq+3*Tgts.sigs'))) * optInputs.penalty;
                    end
                end
                
            elseif optInputs.optType == 1
                for h = 1:length(optInputs.PerTgt)
                    % Sort the lnSa values at each period and calculate the
                    % normal CDF from this 
                    sortedlnSa = [min(sampleSmall(:,h)); sort(sampleSmall(:,h))];
                    norm_cdf = normcdf(sortedlnSa,Tgts.meanReq(h),Tgts.sigs(h));
                    
                    % Calculate the Dn value
                    Dn(h) = max(abs(emp_cdf'-norm_cdf));
                end
                devTotal(j) = sum(Dn);   
            end
            
            
            if (scaleFac(j) > optInputs.maxScale)
                devTotal(j) = devTotal(j) + 1000000;
            end
            
            % Should cause improvement and record should not be repeated
            if (any(optInputs.recID == j))
                devTotal(j) = 100000;
            end
            sampleSmall = sampleSmall(1:end-1,:);
        end
        
        [minDevFinal, minID] = min(devTotal);
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
    notT1 = find(optInputs.PerTgt ~= optInputs.T1);
    sigs = std(sampleSmall);
    meanErr = max(abs(exp(mean(sampleSmall))-Tgts.means)./Tgts.means)*100;
    sigErr = max(abs(sigs(notT1) - Tgts.sigs(notT1))./Tgts.sigs(notT1))*100;
    
    % Display the results
    display(meanErr);
    display(sigErr);
    if meanErr < optInputs.tol && sigErr < optInputs.tol
        display('The percent errors between chosen and target spectra are now within the required tolerances.');
        break;
    end
end

display('100% done');

% Output information
finalRecords = optInputs.recID;
finalScaleFactors = finalScaleFac';




end

