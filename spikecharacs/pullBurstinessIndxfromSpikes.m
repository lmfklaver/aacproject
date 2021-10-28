% Script to pull all the AAC burstiness out of CellExplorer


allBurstRoyer =[];allBurstMizuseki=[]; allBurstDoublets =[];
allSess = [];   allAACIndx = [];
%%
for iSess = [1,2,4,5,15,16,17]%:length(dirN)
    cd(dirN{iSess})
    basepath = cd; basename = bz_BasenameFromBasepath(cd);
    load([basename '.cell_metrics.cellinfo.mat'])
%     spikes = bz_LoadPhy;
    [~,~,aacs] = splitCellTypes(basepath);

    burstRoyer = cell_metrics.burstIndex_Royer2012;
    burstMizuseki = cell_metrics.burstIndex_Mizuseki2012;
    burstDoublets = cell_metrics.burstIndex_Doublets;
   
     
    allBurstRoyer = [allBurstRoyer , burstRoyer(aacs)];
    allBurstMizuseki = [allBurstMizuseki, burstMizuseki(aacs)];
    allBurstDoublets = [allBurstDoublets, burstDoublets(aacs)];
   
    allSess     = [allSess, repmat(iSess,1,length(aacs))];
    allAACIndx  = [allAACIndx, aacs];
end

a = [allSess', allAACIndx',allBurstRoyer',allBurstMizuseki',allBurstDoublets'];