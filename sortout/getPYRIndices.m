% D:\Data\Axoaxo\mouse1\mouse1_180414
%%
PyrIND_concat = [];

for iSess = 1:length(sessions)
    fprintf(num2str(iSess))
    pause
    selecSession = sessions{iSess};
    
      % Load Spikes
    load([folderAACchr filesep selecSession '.cell_metrics.cellinfo.mat']);
    
    for iUnit = 1:length(cell_metrics.putativeCellType)
        PyrIND(iUnit) = strcmp(cell_metrics.putativeCellType{iUnit},{'Pyramidal Cell'});
    end
    PyrIND_concat = [PyrIND_concat, PyrIND];
    
    clear PyrIND
end


%%
PyrIND_concat = logical(PyrIND_concat)