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
% origPerKnown = perKnown;
% if ~any(perKnown == optInputs.T1)
%     perKnown = [perKnown(perKnown<T1) optInputs.T1 perKnown(perKnown>T1)];
% end

perKnown= perKnown(perKnown<=10); % omit long periods for which there is no GMPE

% sigmaKnown = zeros(1,length(perKnown));
% for i = 1:length(perKnown)
%     [tmp sigmaKnown(1,i)] = CB_2008_nga (M_bar, perKnown(i), Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb);
% end

% corrReq = zeros(length(perKnown));
% for i=1:length(perKnown)
%     for j=1:length(perKnown)
% 
%         Ta = perKnown(i);
%         Tb = perKnown(j);
%         
%         if optInputs.cond == 1
%             rec = find(perKnown == optInputs.T1);
%             var1 = sigmaKnown(i)^2;
%             var2 = sigmaKnown(j)^2;
%             varT = sigmaKnown(rec)^2;
% 
%             sigma11 = [var1 baker_jayaram_correlation(Ta, Tb)*sqrt(var1*var2);baker_jayaram_correlation(Ta, Tb)*sqrt(var1*var2) var2];
%             sigma22 = varT;
%             sigma12 = [baker_jayaram_correlation(Ta, optInputs.T1)*sqrt(var1*varT);baker_jayaram_correlation(optInputs.T1, Tb)*sqrt(var2*varT)];
% 
%             sigmaCond = sigma11 - sigma12*inv(sigma22)*(sigma12)';
%             corrReq(i,j) = sigmaCond(1,2)/sqrt(sigmaCond(1,1)*sigmaCond(2,2));
%             
%         elseif optInputs.cond == 0
%             corrReq(i,j) = baker_jayaram_correlation(Ta, Tb);
%         end
%     end
% end

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
    contour(perKnown, perKnown, corrReqSamp);
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
    contour(perKnown, perKnown, corrReq);
    set(gca,'yscale','log','xscale','log'); 
    axis square;
    xlabel('T_1');
    ylabel('T_2');
    title('Baker and Jayaram (2008) conditional correlation contour');
    xlabel('T_1')
    ylabel('T_2')
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    % Error
    diffCorr = corrReqSamp-corrReq;
    figure
    contour(perKnown,perKnown,diffCorr);
    set(gca, 'yscale', 'log','xscale','log');
    axis square;
    title('Difference in the correlation (sample-model)');
    xlabel('T_1 (s)');
    ylabel('T_2 (s)');
    colorbar('YLim',[-1 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
end

