% This code is used to compare the covariance structure of the selected
% ground motions with the covariance structure provided by Baker and
% Jayaram (2008).
%
% Nirmal Jayaram, Ting Lin, Jack W. Baker
% Department of Civil and Environmental Engineering
% Stanford University
% Last Updated: 11 March 2010
%
% Reference manuscripts:
%
% J. W. Baker and Jayaram, N. (2008). Correlation of spectral acceleration 
% values from NGA ground motion models, Earthquake Spectra, 24 (1), 299-317

%% Observed correlations
sampleUse = [];
sampleUse = log(SaKnown(finalRecords,:).*repmat(finalScaleFactors,1,size(SaKnown,2)));
sampleUse = [sampleUse(:,perKnown<optInputs.T1) interp1(perKnown,sampleUse',optInputs.T1)' sampleUse(:,perKnown>optInputs.T1)];
corrReqSamp = zeros(length(perKnownCorr));
for i=1:length(perKnownCorr)
    for j=1:length(perKnownCorr)
        corrMatrix = corrcoef((sampleUse(:,i)),(sampleUse(:,j)));
        corrReqSamp(i,j) = corrMatrix(1,2);
    end
end
   
if showPlots == 1
%% contour plots

    figure
    contour(perKnownCorr, perKnownCorr, corrReqSamp);
    set(gca,'yscale','log','xscale','log'); 
    axis square;
    xlabel('T_1');
    ylabel('T_2');
    title('Sample correlation coefficients contour');
    xlabel('T_1 (s)');
    ylabel('T_2 (s)');
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)

    figure
    contour(perKnownCorr, perKnownCorr, corrReq);
    set(gca,'yscale','log','xscale','log'); 
    axis square;
    xlabel('T_1 (s)');
    ylabel('T_2 (s)');
    title('Target correlation coefficients contour');
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    % Error
    diffCorr = corrReqSamp-corrReq;
    figure
    contour(perKnownCorr,perKnownCorr,diffCorr);
    set(gca, 'yscale', 'log','xscale','log');
    axis square;
    title('Difference in the correlation (sample-target)');
    xlabel('T_1 (s)');
    ylabel('T_2 (s)');
    colorbar('YLim',[-1 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
end

