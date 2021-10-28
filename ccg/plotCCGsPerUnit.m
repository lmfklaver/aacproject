%% CCGs for ccginout
%  Get the CCG of the spikes outside the pulse (NB , CCG function in
%  buzcode, conflicts with CCG in CellExplorer, so remove CellExplorer
%  from path)


% Make figure and save basename.ccg.mat
saveFig    = true;
doPlot     = true;

savePath = cd;
saveName = [basename '_CCGs_' date];

%%
% plot CCG
ccg = ccginout.ccgOUT;
t = ccginout.t;

if doPlot
    
    for iUnit = 1:size(ccg,2)
        numPanels = ceil(sqrt(size(ccg,2)));
        
        figure
        set(gcf,'Position',[50 50 1200 800]);
        set(gcf,'PaperOrientation','landscape');
        
        for nPlot = 1:size(ccg,3)
            subplot(numPanels,numPanels,nPlot)
            if iUnit == nPlot
                bar(t,ccg(:,iUnit,nPlot),'k')
            else
                bar(t,ccg(:,iUnit,nPlot))
            end
            title(num2str(spikes.cluID(nPlot)));
        end
        
        %             s1 = suptitle([num2str(iSess) '_' num2str(spikes.UID(iUnit))]);
        %             set(s1, 'Interpreter','none')
        
        if saveFig
            unitStr = num2str(spikes.UID(iUnit));
            print(gcf,[unitStr '.pdf'],'-dpdf','-bestfit')
            append_pdfs(fullfile(savePath,[saveName '.pdf']),[unitStr '.pdf'])
            delete([unitStr '.pdf'])
            close
        end
    end
end


