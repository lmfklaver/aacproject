
%gives phase (x) by power bin (y) by rate (z)?


filteredLFP = bz_Filter(lfp,'filter','butter','passband',[1 100],'order', 3);
[PowerPhaseRatemap,spikebinIDs] = bz_PowerPhaseRatemap(spikes,filteredLFP);

subplot(2,2,1)
    imagesc(PowerPhaseRatemap.phasebins,PowerPhaseRatemap.powerbins,...
        PowerPhaseRatemap.meanrate)
    hold on
    imagesc(PowerPhaseRatemap.phasebins+2*pi,PowerPhaseRatemap.powerbins,...
        PowerPhaseRatemap.meanrate)
    plot(linspace(-pi,3*pi,100),cos(linspace(-pi,3*pi,100)),'k')
    xlim([-pi 3*pi])
    axis xy
    colorbar  
    xlabel('Phase');
    ylabel('Power (Z(log))')
    
subplot(2,2,3)
    bar(PowerPhaseRatemap.powerbins,PowerPhaseRatemap.powerdist)
    xlabel('Power (Z(log))');ylabel('Time (s)')
    box off
    axis tight
    
    