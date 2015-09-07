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

%% Estimate conditional covariances from the Baker and Jayaram (2008)
% model

% Modify perTgt to include T1
if ~any(perKnown == optInputs.T1)
    perKnown1 = [perKnown(perKnown<T1) optInputs.T1 perKnown(perKnown>T1)];
else
    perKnown1 = perKnown;
end

perKnown1= perKnown1(perKnown1<=10); % omit long periods for which there is no GMPE

sigmaKnown = zeros(1,length(perKnown1));
for i = 1:length(perKnown1)
    [tmp sigmaKnown(1,i)] = CB_2008_nga (M_bar, perKnown1(i), Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb);
end


corrReq = zeros(length(perKnown1));
corrReqSamp = zeros(length(perKnown1));
for i=1:length(perKnown1)
    for j=1:length(perKnown1)
        
        Ta = perKnown1(i);
        Tb = perKnown1(j);
        
        rec = find(perKnown1 == optInputs.T1);
        var1 = sigmaKnown(i)^2;
        var2 = sigmaKnown(j)^2;
        varT = sigmaKnown(rec)^2;
        
        
        sigma11 = [var1 baker_jayaram_correlation(Ta, Tb)*sqrt(var1*var2);baker_jayaram_correlation(Ta, Tb)*sqrt(var1*var2) var2];
        sigma22 = varT;
        sigma12 = [baker_jayaram_correlation(Ta, optInputs.T1)*sqrt(var1*varT);baker_jayaram_correlation(optInputs.T1, Tb)*sqrt(var2*varT)];
        
        sigmaCond = sigma11 - sigma12*inv(sigma22)*(sigma12)';
        
        corrReq(i,j) = sigmaCond(1,2)/sqrt(sigmaCond(1,1)*sigmaCond(2,2));
        
    end
end

if showPlots == 1
    figure
    imagesc(perKnown1,perKnown1,corrReq)
    title('Baker and Jayaram (2008) conditional correlations');
    xlabel('T_1 (s)');
    ylabel('T_2 (s)');
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    set(gca, 'YDir', 'normal')
end

%% Observed correlations

sampleUse = [];
sampleUse = log(SaKnown(finalRecords,:).*repmat(finalScaleFactors,1,size(SaKnown,2)));
sampleUse = [sampleUse(:,perKnown<optInputs.T1) interp1(perKnown,sampleUse',optInputs.T1)' sampleUse(:,perKnown>optInputs.T1)];

for i=1:length(perKnown1)
    for j=1:length(perKnown1)
        corrMatrix = corrcoef((sampleUse(:,i)),(sampleUse(:,j)));
        corrReqSamp(i,j) = corrMatrix(1,2);
    end
end
   
if showPlots == 1
    figure
    imagesc(perKnown1,perKnown1,corrReqSamp)
    title('Sample correlations');
    xlabel('T_1 (s)');
    ylabel('T_2 (s)');
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    set(gca, 'YDir', 'normal')
end

%% Error
diffCorr = corrReqSamp-corrReq;
% Display minimum and maximum difference in correlations
% display(min(abs(diffCorr(:))));
% display(max(abs(diffCorr(:))));
if showPlots == 1
    figure
    imagesc(perKnown1,perKnown1,diffCorr)
    title('Difference in the correlation (sample-model)');
    xlabel('T_1 (s)');
    ylabel('T_2 (s)');
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    set(gca, 'YDir', 'normal')
end

%% contour plot

if showPlots == 1
    figure
    contour(perKnown1, perKnown1, corrReqSamp);
    set(gca,'yscale','log','xscale','log'); 
    axis square;
    xlabel('T_1');
    ylabel('T_2');
    title('Sample correlation contour');
    xlabel('T_1')
    ylabel('T_2')
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)

    figure
    contour(perKnown1, perKnown1, corrReq);
    set(gca,'yscale','log','xscale','log'); 
    axis square;
    xlabel('T_1');
    ylabel('T_2');
    title('Baker and Jayaram (2008) conditional correlation contour');
    xlabel('T_1')
    ylabel('T_2')
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
end

