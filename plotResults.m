%% Produce a number of plots of target and selected response spectra 

    % Variables used here
    
    % SaKnown    : As before, it contains the response spectra of all the
    %              available ground motions (N*P matrix) - N ground motions,
    %              P periods
    % simulatedSpectra  : a matrix of simulated response spectra defined
    %              only at PerTgt
    % sampleBig  : Same as SaKnown, but is only defined at PerTgt, the
    %              periods at which the target response spectrum properties
    %              are computed
    % sampleSmall: The response spectra of the selected ground motions,
    %              defined at PerTgt
    % meanReq    : Target mean for the (log) response spectrum
    % covReq     : Target covariance for the (log) response spectrum

    
    % Plot simulated response spectra  
    figure
    loglog(optInputs.PerTgt, exp(Tgts.meanReq), '-r', 'linewidth', 3)
    hold on
    loglog(optInputs.PerTgt, exp(Tgts.meanReq + 1.96*sqrt(diag(Tgts.covReq))'), '--r', 'linewidth', 3)
    loglog(optInputs.PerTgt, simulatedSpectra','k');
    loglog(optInputs.PerTgt, exp(Tgts.meanReq - 1.96*sqrt(diag(Tgts.covReq))'), '--r', 'linewidth', 3)
    axis([min(optInputs.PerTgt) max(optInputs.PerTgt) 1e-2 5])
    xlabel('T (s)')
    ylabel('S_a (g)')
    legend('Median response spectrum','2.5 and 97.5 percentile response spectra','Response spectra of simulated ground motions')
    title('Response spectra of simulated ground motions')
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    % Plot selected response spectra
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

    % Target, initial, and finally selected medians
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
    
    % Target, initial, and finally selected standard deviations
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


    % Compute sample correlations from selected spectra
    sampleUse = log(SaKnown(finalRecords,:).*repmat(finalScaleFactors,1,size(SaKnown,2)));
    sampleUse = [sampleUse(:,perKnown<optInputs.T1) interp1(perKnown,sampleUse',optInputs.T1)' sampleUse(:,perKnown>optInputs.T1)];
    corrReqSamp = zeros(length(perKnownCorr));
    for i=1:length(perKnownCorr)
        for j=1:length(perKnownCorr)
            corrMatrix = corrcoef((sampleUse(:,i)),(sampleUse(:,j)));
            corrReqSamp(i,j) = corrMatrix(1,2);
        end
    end
    
    % Contours of correlations from selected spectra
    figure
    contour(perKnownCorr, perKnownCorr, corrReqSamp);
    set(gca,'yscale','log','xscale','log');
    axis square;
    xlabel('T_1');
    ylabel('T_2');
    title('Sample correlation coefficients contour');
    xlabel('T_1 (s)')
    ylabel('T_2 (s)')
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    % Contours of target correlations
    figure
    contour(perKnownCorr, perKnownCorr, corrReq);
    set(gca,'yscale','log','xscale','log');
    axis square;
    xlabel('T_1 (s)');
    ylabel('T_2 (s)');
    title('Target correlation coefficients contour');
    colorbar('YLim',[0 1]);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    % Difference between target and sample correlations
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
    

