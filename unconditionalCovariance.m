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

%% Estimate unconditional covariances from the Baker and Jayaram (2008)
% model

corrReq = zeros(length(perKnown));
corrReqSamp = zeros(length(perKnown));
for i=1:length(perKnown)
    for j=1:length(perKnown)
        
        T1 = perKnown(i);
        T2 = perKnown(j);
        corrReq(i,j) = baker_jayaram_correlation(T1, T2);
        
    end
end

%% Observed correlations

sampleUse = [];
sampleUse = [sampleUse;log(SaKnown(finalRecords,:).*repmat(finalScaleFactors,1,size(SaKnown,2)))];

for i=1:length(perKnown)
    for j=1:length(perKnown)
        corrMatrix = corrcoef((sampleUse(:,i)),(sampleUse(:,j)));
        corrReqSamp(i,j) = corrMatrix(1,2);
    end
end

%% contour plots

if showPlots == 1
    figure
    contour(perKnown, perKnown, corrReq);
    set(gca,'yscale','log','xscale','log'); 
    xlabel('T_1');
    ylabel('T_2');
    title('Baker and Jayaram (2008) correlation contour');
    xlabel('T_1')
    ylabel('T_2')
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    figure
    contour(perKnown, perKnown, corrReqSamp);
    set(gca,'yscale','log','xscale','log'); 
    xlabel('T_1');
    ylabel('T_2');
    title('Sample correlation contour');
    xlabel('T_1')
    ylabel('T_2')
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18) 
    
    % Error
    diffCorr = corrReqSamp-corrReq;
    figure
    imagesc(perKnown,perKnown,diffCorr)
    title('Difference in the correlation (sample-model)');
    xlabel('T_1 (s)');
    ylabel('T_2 (s)');
    colorbar('YLim', [-1 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    set(gca,'yscale','log','xscale','log'); 
end
