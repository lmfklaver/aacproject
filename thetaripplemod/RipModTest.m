for iSess = 1:length(dirN)

    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(basepath);
    
    load([basename '.spikes.cellinfo.mat'])
    ripmod = getRipMod_CCG(basepath);
    
end

    
    
   
    