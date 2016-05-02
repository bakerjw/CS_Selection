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
%            indTcond   : The index of Tcond, the conditioning period
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

if selectionParams.cond == 0 && selectionParams.isScaled
    display('The algorithm is slower when scaling is used');
end
if selectionParams.optType == 1
    display('The algorithm is slower when optimizing with the KS-test Dn statistic');
end


% Initialize scale factor vectors if possible

if selectionParams.isScaled == 0 % no scaling so set scale factors = 1
    scaleFac = ones(selectionParams.nBig,1);
    IMs.scaleFac = ones(selectionParams.nGM,1);
elseif selectionParams.isScaled && selectionParams.cond % Sa(Tcond) scaling
    scaleFac = exp(selectionParams.lnSa1)./exp(IMs.sampleBig(:,selectionParams.indTcond));
    % get indices of ground motions with allowable scale factors, for further consideration
    idxAllow = find(scaleFac < selectionParams.maxScale);
end

hw = waitbar(0,'Optimizing ground motion selection');

for k=1:selectionParams.nLoop % Number of passes
    for i=1:selectionParams.nGM % consider replacing each ground motion in the selected set

        sampleSmall(i,:) = []; % remove initially selected record to be replaced
        IMs.recID(i,:) = []; 
        
        % if scaling with unconditional selection, compute scale factors
        if selectionParams.isScaled && selectionParams.cond == 0
            scaleFac = compute_scale_factor(IMs.sampleBig, sampleSmall, targetSa.meanReq, targetSa.stdevs, selectionParams.weights, selectionParams.maxScale);
            % get indices of ground motions with allowable scale factors, for further consideration
            idxAllow = find(scaleFac < selectionParams.maxScale);
        end
        
        % Try to add a new spectrum to the subset list
        devTotal = 1000000 * ones(selectionParams.nBig,1); % initialize to large errors, and recompute for allowable records
        for j = 1:length(idxAllow) 
            
            if ~any(IMs.recID == idxAllow(j)) % if this candidate is not already in the set 
                testSpectra = [sampleSmall; IMs.sampleBig(idxAllow(j),:)+log(scaleFac(idxAllow(j)))]; % add candidate to set
                devTotal(idxAllow(j)) = compute_spectrum_error(selectionParams, targetSa, testSpectra);
            end
            
        end
        
        [~ , minID] = min(devTotal);
        % Add new element in the right slot
        IMs.scaleFac(i) = scaleFac(minID);
        IMs.recID = [IMs.recID(1:i-1); minID; IMs.recID(i:end)];
        sampleSmall = [sampleSmall(1:i-1,:); IMs.sampleBig(minID,:)+log(scaleFac(minID)); sampleSmall(i:end,:)];
        
        waitbar(((k-1)*selectionParams.nGM + i)/(selectionParams.nLoop*selectionParams.nGM)); % update waitbar
    end
    
    % check whether results are within tolerance, and stop optimization if so
    if within_tolerance(sampleSmall, targetSa, selectionParams)
        break;
    end
end

close(hw); % close waitbar

% Save final selection for output
IMs.sampleSmall = sampleSmall;


end

