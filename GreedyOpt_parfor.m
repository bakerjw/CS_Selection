function [ sampleSmall, finalRecords, finalScaleFactors ] = GreedyOpt_parfor( optInputs, Tgts, IMs )
% Greedy optimization with parallelization

% define all variables within structures to eliminate unnecessary
% communication overhead (try and reduce information sent to each worker?)
sampleSmall = IMs.sampleSmall;
sampleBig   = IMs.sampleBig;
isScaled    = optInputs.isScaled;
cond        = optInputs.cond;
PerTgt      = optInputs.PerTgt;
T1          = optInputs.T1;

if cond == 1
    lnSa1       = optInputs.lnSa1;
else
    lnSa1       = [];
end


display('Please wait...This algorithm takes a few minutes depending on the number of records to be selected');
if optInputs.cond == 0
    display('The algorithm is slower when scaling is used');
end

% open parallel pool workers
parobj = parpool(3);

% begin optimization
% for k=1:optInputs.nLoop % Number of passes
    
    for i=1:optInputs.nGM % Selects nSelect ground motions
        display([num2str(round(((k-1)*optInputs.nGM + i-1)/(optInputs.nLoop*optInputs.nGM)*100)) '% done']);
        
        minDev = zeros(optInputs.nBig,1);
        sampleSmall1 = sampleSmall([1:i-1, i+1:end],:);
        optInputs.recID(i,:) = [];
        
        scaleFac = zeros(optInputs.nBig,1);
        
        % Try to add a new spectra to the subset list on multiple workers
        parfor j = 1:optInputs.nBig
            if isScaled == 1
                if cond == 1
                    scaleFac(j) = exp(lnSa1)/exp(sampleBig(j,PerTgt == T1));
                elseif cond == 0
                    [scaleFac(j), devTotal] = bestScaleFactor(sampleBig(j,:),sampleSmall1,Tgts.meanReq,sqrt(diag(Tgts.covReq))',optInputs.weights,optInputs.maxScale);
                end
                sampleSmallTemp = [sampleSmall1;sampleBig(j,:)+log(scaleFac(j))];
            else
                sampleSmallTemp = [sampleSmall1;sampleBig(j,:)];
                scaleFac(j) = 1;
            end
            
            if cond == 1 | (cond == 0 && isScaled == 0)
                % Compute deviations from target
                devMean = mean(sampleSmallTemp) - Tgts.meanReq;
                devSig = std(sampleSmallTemp) - sqrt(diag(Tgts.covReq))';
                devTotal = optInputs.weights(1) * sum(devMean.^2) + optInputs.weights(2) * sum(devSig.^2);
            end
            
            % Penalize bad spectra (set penalty to zero if this is not required)
            if optInputs.penalty ~= 0
                for m=1:size(sampleSmallTemp,1)
                    devTotal = devTotal + sum(abs(exp(sampleSmallTemp(m,:))>exp(Tgts.meanReq+3*sqrt(diag(Tgts.covReq))'))) * optInputs.penalty;
                end
            end
            
            if (scaleFac(j) > optInputs.maxScale)
                devTotal = devTotal + 1000000;
            end
            
%             % Should cause improvement and record should not be repeated
%             if (devTotal < minDev && ~any(optInputs.recID == j))
%                 minID = j;
%                 minDev = devTotal;
%             end
             if (~any(optInputs.recID == j))
                 minDev(j) = min(devTotal);
             else 
                 minDev(j) = 100000;
             end
            sampleSmallNew = sampleSmallTemp(1:end-1,:);
        end
        
        [minDevFinal, minID] = min(minDev);
       
        % Add new element in the right slot
        if optInputs.isScaled == 1
            finalScaleFac(i) = scaleFac(minID);
        else
            finalScaleFac(i) = 1;
        end
        sampleSmall = [sampleSmall(1:i-1,:);IMs.sampleBig(minID,:)+log(scaleFac(minID));sampleSmall(i:end,:)];
        optInputs.recID = [optInputs.recID(1:i-1);minID;optInputs.recID(i:end)];
        
    end
    
    % can the optimization be stopped after one loop?
    percentMean = max(abs(exp(mean(sampleSmall))-exp(Tgts.meanReq))./exp(Tgts.meanReq))*100;
    
    if optInputs.cond == 1
        T1_index = find(optInputs.PerTgt == optInputs.T1);
        sigs = sqrt(diag(Tgts.covReq))';
        sampleSigs = std(sampleSmall);
        percentSig = max(abs(sampleSigs([1:T1_index-1,T1_index+1:end])-sigs([1:T1_index-1,T1_index+1:end]))./sigs([1:T1_index-1,T1_index+1:end]))*100;
    else
        percentSig = max(abs(std(sampleSmall)-sqrt(diag(Tgts.covReq))')./sqrt(diag(Tgts.covReq))')*100;
    end
    
    if percentMean < optInputs.tol && percentSig < optInputs.tol
        display('The percent errors between chosen and target spectra are now within the required tolerances.');
        break
    end
% end

delete(parobj);

display('100% done');

% Output information
finalRecords = optInputs.recID;
finalScaleFactors = finalScaleFac';

end

