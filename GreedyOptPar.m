function [ sampleSmall, finalRecords, finalScaleFactors ] = GreedyOptPar( optInputs, Tgts, IMs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
% period) and initialize a vector of Dn values
emp_cdf = 0;
if optInputs.optType == 1
    emp_cdf = linspace(0,1,optInputs.nGM+1);
end

numWorkers = 2;
parobj = parpool(numWorkers);

% Initialize scale factor vector
scaleFac = ones(optInputs.nBig,1);
for k=1:optInputs.nLoop % Number of passes
    
    for i=1:optInputs.nGM % Selects nGM ground motions
        display([num2str(round(((k-1)*optInputs.nGM + i-1)/(optInputs.nLoop*optInputs.nGM)*100)) '% done']);
        
        devTotal = zeros(optInputs.nBig,1);
        sampleSmall(i,:) = [];
        optInputs.recID(i,:) = [];
        
        if optInputs.isScaled == 1
            if optInputs.cond == 1
                scaleFac = exp(optInputs.lnSa1)./exp(IMs.sampleBig(:,optInputs.rec));
            elseif optInputs.cond == 0
                [scaleFac, devTotal] = bestScaleFactor(IMs.sampleBig, sampleSmall, Tgts.meanReq, Tgts.stdevs, optInputs.weights, optInputs.maxScale);
            end
        end
        
        % Try to add a new spectra to the subset list
        % new function 
        [devTotal] = ParLoop(devTotal, scaleFac, optInputs, sampleSmall, IMs.sampleBig, Tgts.meanReq,...
                             Tgts.stdevs, emp_cdf);                                
                                
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
    if optInputs.optType == 0
        notT1 = find(optInputs.PerTgt ~= optInputs.PerTgt(optInputs.rec));
        stdevs = std(sampleSmall);
        meanErr = max(abs(exp(mean(sampleSmall))-Tgts.means)./Tgts.means)*100;
        stdErr = max(abs(stdevs(notT1) - Tgts.stdevs(notT1))./Tgts.stdevs(notT1))*100;
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
finalScaleFactors = finalScaleFac';

delete(parobj);
end

