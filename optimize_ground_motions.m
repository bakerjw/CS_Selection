function [ IMs ] = optimize_ground_motions( selectionParams, targetSa, IMs )
% This function will perform a greedy optimization on a set of ground
% motions using the sum of squared errors approach to check the set of
% selected ground motions against target means and variances
%
% selectionParams: input values needed to run the optimization function
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
% targetSa    :  The target values (means and covariances) being matched
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

% sampleSmall changes size throughout the optimization. Redfine sampleSmall
% here. sampleSmall is returned as a new variable, not within IMs
sampleSmall = IMs.sampleSmall;

display('Please wait...This algorithm takes a few minutes depending on the number of records to be selected');
if selectionParams.cond == 0
    display('The algorithm is slower when scaling is used');
end
if selectionParams.optType == 1
    display('The algorithm is slower when optimizing with the KS-test Dn statistic');
end

% if optimizing with KS-test statistic, first calculate the emperical CDF
% values (which are the same at each period)
if selectionParams.optType == 1
    emp_cdf = linspace(0,1,selectionParams.nGM+1);
end

% Initialize scale factor vector for each iteration and a scale factor
% vector for the final selected records 
scaleFac = ones(selectionParams.nBig,1);
IMs.scaleFac = ones(selectionParams.nGM,1);

hw = waitbar(0,'Optimizing ground motion selection');

for k=1:selectionParams.nLoop % Number of passes
    for i=1:selectionParams.nGM % Selects nGM ground motions
        waitbar(((k-1)*selectionParams.nGM + i-1)/(selectionParams.nLoop*selectionParams.nGM)); % update percentage completion 

        devTotal = zeros(selectionParams.nBig,1); % initialize measure of deviation for each selected ground motion
        sampleSmall(i,:) = []; % remove initially selected record to be replaced
        IMs.recID(i,:) = []; 
        
        if selectionParams.isScaled == 1
            if selectionParams.cond == 1
                scaleFac = exp(selectionParams.lnSa1)./exp(IMs.sampleBig(:,selectionParams.indT1));
            elseif selectionParams.cond == 0
                [scaleFac, devTotal] = bestScaleFactor(IMs.sampleBig, sampleSmall, targetSa.meanReq, targetSa.stdevs, selectionParams.weights, selectionParams.maxScale);
            end
        end
        
        % Try to add a new spectrum to the subset list
        for j = 1:selectionParams.nBig
            sampleSmall = [sampleSmall;IMs.sampleBig(j,:)+log(scaleFac(j))];
            % Calculate the candidate's deviation and store in devTotal 
            if selectionParams.optType == 0
                if selectionParams.cond == 1 || (selectionParams.cond == 0 && selectionParams.isScaled == 0)
                    % Compute deviations from target
                    devMean = mean(sampleSmall) - targetSa.meanReq;
                    devSig = std(sampleSmall) - targetSa.stdevs;
                    devTotal(j) = selectionParams.weights(1) * sum(devMean.^2) + selectionParams.weights(2) * sum(devSig.^2);
                end
                % Penalize bad spectra (set penalty to zero if this is not required)
                if selectionParams.penalty ~= 0
                    for m=1:size(sampleSmall,1)
                        devTotal(j) = devTotal(j) + sum(abs(exp(sampleSmall(m,:))>exp(targetSa.meanReq+3*targetSa.stdevs'))) * selectionParams.penalty;
                    end
                end
                
            elseif selectionParams.optType == 1   
                sortedlnSa = [min(sampleSmall); sort(sampleSmall)];
                norm_cdf = normcdf(sortedlnSa,repmat(targetSa.meanReq,selectionParams.nGM+1,1),repmat(targetSa.stdevs,selectionParams.nGM+1,1));
                Dn = max(abs(repmat(emp_cdf',1,length(selectionParams.TgtPer)) - norm_cdf));
                devTotal(j) = sum(Dn);   
            end
            
            % Scale factor should not exceed the maximum
            if (scaleFac(j) > selectionParams.maxScale)
                devTotal(j) = 1000000;
            end
            
            % Record should not be repeated
            if (any(IMs.recID == j))
                devTotal(j) = 1000000;
            end
            sampleSmall = sampleSmall(1:end-1,:);
        end
        
        [~ , minID] = min(devTotal);
        % Add new element in the right slot
        if selectionParams.isScaled == 1
            IMs.scaleFac(i) = scaleFac(minID);
        else
            IMs.scaleFac(i) = 1;
        end
        sampleSmall = [sampleSmall(1:i-1,:);IMs.sampleBig(minID,:)+log(scaleFac(minID));sampleSmall(i:end,:)];
        IMs.recID = [IMs.recID(1:i-1);minID;IMs.recID(i:end)];
    end
    
    % Can the optimization be stopped after this loop based on the user
    % specified tolerance? Recalculate new means and standard deviations of
    % new sampleSmall and then recalculate new maximum percent errors of
    % means and standard deviations
    if selectionParams.optType == 0
        stdevs = std(sampleSmall);
        meanErr = max(abs(exp(mean(sampleSmall))- exp(targetSa.meanReq))./exp(targetSa.meanReq))*100;
        stdErr = max(abs(stdevs(1:end ~= selectionParams.indT1) - targetSa.stdevs(1:end ~= selectionParams.indT1))./targetSa.stdevs(1:end ~= selectionParams.indT1))*100;
        fprintf('Max (across periods) error in median = %3.1f percent \n', meanErr); 
        fprintf('Max (across periods) error in standard deviation = %3.1f percent \n \n', stdErr);
        
        % If error is now within the tolerance, break out of the
        % optimization
        if meanErr < selectionParams.tol && stdErr < selectionParams.tol
            display('The percent errors between chosen and target spectra are now within the required tolerances.');
            break;
        end
    end
    
    fprintf('End of loop %1.0f of %1.0f \n', k, selectionParams.nLoop) 
end

close(hw); % close waitbar

% Save final selection for output
IMs.sampleSmall = sampleSmall;


end

