%% inclusion for PYRs
cd(basepath)
basename = bz_BasenameFromBasepath(basepath)

[clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityKilosort([basepath]);% filesep 'Kilosort2']);
save([basename '.sortinquality.cellinfo.mat'],'clusterIDs','unitQuality','contaminationRate')
cd(basepath)
load([basename '.cell_metrics.cellinfo.mat'])
load([basename '.spikes.cellinfo.mat'])


%%
idx = ismember(clusterIDs, cell_metrics.cluID+1);
contaR = contaminationRate(idx);
unQ = unitQuality(idx);
%%
figure

for iUnit= 1:spikes.numcells
    if contains(cell_metrics.putativeCellType{iUnit},'Interneuron')
        if cell_metrics.firingRate(iUnit)>10
            plot(contaR(iUnit),unQ(iUnit),'bo','MarkerFaceColor','b')
        else
            plot(contaR(iUnit),unQ(iUnit),'bo')
        end
        hold on
        
    else
        if cell_metrics.firingRate(iUnit)>10
            plot(contaR(iUnit),unQ(iUnit),'r^','MarkerFaceColor','r')
        else
            plot(contaR(iUnit),unQ(iUnit),'r^')
        end
        
        hold on
    end
end

box off 
xlabel('contamination Rate')
ylabel('unit Quality')
title([basename],'Interpreter','none')

savefig([basename '_sortingQuality'])