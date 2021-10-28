%% Get PETH of All Cells in One plot

% This is to get an overview plot of all PETHs per session in one plot
    
launchDirNforAACSessions

%%
for iSess = 17
    
    % Go to folder and load in spikes and optostim
    
    cd(dirN{iSess})
    
    basepath = cd; 
    basename = bz_BasenameFromBasepath(cd);
    
    spikes = bz_LoadPhy;
    
    load([basename '.optoStim.manipulation.mat'])
    pulseEpochs = optoStim.timestamps;
    
    % Set parameters for PETHs
   [peth] = getPETH_epochs(basepath,'epochs',pulseEpochs);
     
    ratePulse   = peth.rate;
    countPulse  = peth.count;
    timeEdges   = peth.timeEdges; 
    
        
    %%  plot PETH
    figure
    
    numPanels = ceil(sqrt(length(spikes.times)));
    set(gcf,'PaperOrientation','Landscape')
    
    for iUnit = 1:length(spikes.times)
        subplot(numPanels,numPanels,iUnit)
        histogram('BinEdges',timeEdges, 'BinCounts',ratePulse(iUnit,:))
        hold on
        box off
        set(gca,'TickDir','out')
        title({['Unit ' num2str(iUnit)]})
        xlabel('time(s)')
        ylabel('spikes/s')
%         xlim(timwin)
    end
    
    print(gcf,[num2str(iSess) '_' basename '_PETH_all.pdf'],'-dpdf','-bestfit')
end
