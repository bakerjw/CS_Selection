%% Spectra Plots

    % Variables used here
    
    % SaKnown    : As before, it contains the response spectra of all the
    %              available ground motions (N*P matrix) - N ground motions,
    %              P periods
    % gm         : gm is a matrix of simulated response spectra defined
    %              only at PerTgt
    % sampleBig  : Same as SaKnown, but is only defined at PerTgt, the
    %              periods at which the target response spectrum properties
    %              are computed
    % sampleSmall: The response spectra of the selected ground motions,
    %              defined at PerTgt
    % meanReq    : Target mean for the (log) response spectrum
    % covReq     : Target covariance for the (log) response spectrum

    
    % Plot simulated response spectra -- move with the rest of the figures 
    figure
    loglog(optInputs.PerTgt, exp(Tgts.meanReq), '-r', 'linewidth', 3)
    hold on
    loglog(optInputs.PerTgt, exp(Tgts.meanReq + 1.96*sqrt(diag(Tgts.covReq))'), '--r', 'linewidth', 3)
    loglog(optInputs.PerTgt,gm','k');
    loglog(optInputs.PerTgt, exp(Tgts.meanReq - 1.96*sqrt(diag(Tgts.covReq))'), '--r', 'linewidth', 3)
    axis([min(optInputs.PerTgt) max(optInputs.PerTgt) 1e-2 5])
    xlabel('T (s)')
    ylabel('S_a (g)')
    legend('Median response spectrum','2.5 and 97.5 percentile response spectra','Response spectra of simulated ground motions')
    title('Response spectra of simulated ground motions')
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    % Plot target and simulated means
    figure
    loglog(optInputs.PerTgt,exp(Tgts.meanReq))
    hold on
    loglog(optInputs.PerTgt,exp(mean(log(gm))),'--')
    axis([min(optInputs.PerTgt) max(optInputs.PerTgt) 1e-2 5])
    xlabel('T (s)')
    ylabel('Median S_a (g)')
    legend('exp(Target mean lnS_a)','exp(Mean of simulated lnS_a)')
    title('Target and simulated exponential logarithmic means (i.e., medians)')
    set(findall(gcf,'-property','FontSize'),'FontSize',18)

    % Plot target and simulated standard deviations
    figure
    semilogx(optInputs.PerTgt,sqrt(diag(Tgts.covReq))')
    hold on
    semilogx(optInputs.PerTgt,std(log(gm)),'--')
    axis([min(optInputs.PerTgt) max(optInputs.PerTgt) 0 1])
    xlabel('T (s)')
    ylabel('Standard deviation of lnS_a')
    legend('Target standard deviation of lnS_a','Standard deviation of simulated lnS_a')
    title('Target and simulated logarithmic standard deviations')
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    % Plot at all periods
    figure
    loglog(optInputs.PerTgt, exp(Tgts.meanReq), 'b', 'linewidth', 3)
    hold on
    loglog(optInputs.PerTgt, exp(Tgts.meanReq + 1.96*sqrt(diag(Tgts.covReq))'), '--b', 'linewidth', 3)
    loglog(perKnown,SaKnown(finalRecords,:).*repmat(finalScaleFactors,1,size(SaKnown,2)),'k');
    loglog(optInputs.PerTgt, exp(Tgts.meanReq - 1.96*sqrt(diag(Tgts.covReq))'), '--b', 'linewidth', 3)
    axis([min(optInputs.PerTgt) max(optInputs.PerTgt) 1e-2 5])
    xlabel('T (s)');
    ylabel('S_a (g)');
    legend('Median response spectrum','2.5 and 97.5 percentile response spectra','Response spectra of selected ground motions');
    title ('Response spectra of selected ground motions');
    set(findall(gcf,'-property','FontSize'),'FontSize',18)

    % Sample, original sample, and target means
    figure
    loglog(optInputs.PerTgt,Tgts.means,'k','linewidth',1)
    hold on
    loglog(optInputs.PerTgt, origMeans,'r*', 'linewidth',2)
    loglog(optInputs.PerTgt,exp(mean(IMs.sampleSmall)),'b--','linewidth',1)
    axis([min(optInputs.PerTgt) max(optInputs.PerTgt) 1e-2 5])
    xlabel('T (s)');
    ylabel('Median S_a (g)');
    legend('exp(Target mean lnS_a)','exp(Mean of originally selected lnS_a', 'exp(Mean of selected lnS_a)');
    title('Target and sample exponential logarithmic means (i.e., medians)')
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    % Sample, original sample, and target standard deviations
    figure
    semilogx(optInputs.PerTgt,Tgts.stdevs,'k','linewidth',1)
    hold on
    semilogx(optInputs.PerTgt, origStdevs,'r*','linewidth',2)
    semilogx(optInputs.PerTgt,std(IMs.sampleSmall),'b--','linewidth',1)
    axis([min(optInputs.PerTgt) max(optInputs.PerTgt) 0 1])
    xlabel('T (s)');
    ylabel('Standard deviation of lnS_a');
    legend('Target standard deviation of lnS_a','Standard deviation of originally selected lnS_a','Standard deviation of selected lnS_a');
    title('Target and sample logarithmic standard deviations')
    set(findall(gcf,'-property','FontSize'),'FontSize',18)


%% Correlation Comparison
% This portion of the plot script is used to compare the covariance
% structure of the selected ground motions with the covariance structure
% provided by Baker and Jayaram (2008).
%
% Nirmal Jayaram, Ting Lin, Jack W. Baker Department of Civil and
% Environmental Engineering Stanford University Last Updated: 11 March 2010
%
% Reference manuscripts:
%
% J. W. Baker and Jayaram, N. (2008). Correlation of spectral acceleration
% values from NGA ground motion models, Earthquake Spectra, 24 (1), 299-317

if (checkCorr)
    %% Observed correlations
%     sampleUse = [];
    sampleUse = log(SaKnown(finalRecords,:).*repmat(finalScaleFactors,1,size(SaKnown,2)));
    sampleUse = [sampleUse(:,perKnown<optInputs.T1) interp1(perKnown,sampleUse',optInputs.T1)' sampleUse(:,perKnown>optInputs.T1)];
    corrReqSamp = zeros(length(perKnownCorr));
    for i=1:length(perKnownCorr)
        for j=1:length(perKnownCorr)
            corrMatrix = corrcoef((sampleUse(:,i)),(sampleUse(:,j)));
            corrReqSamp(i,j) = corrMatrix(1,2);
        end
    end
    
    %% contour plots
    
    figure
    contour(perKnownCorr, perKnownCorr, corrReqSamp);
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
    contour(perKnownCorr, perKnownCorr, corrReq);
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
    contour(perKnownCorr,perKnownCorr,diffCorr);
    set(gca, 'yscale', 'log','xscale','log');
    axis square;
    title('Difference in the correlation (sample-model)');
    xlabel('T_1 (s)');
    ylabel('T_2 (s)');
    colorbar('YLim',[-1 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
end

