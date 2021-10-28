% Script to get all the spike widths out for all the units


allSpkWidth =[];  allSess = [];   allAACIndx = [];
%%
for iSess = [1,2,4,5,15,16,17]%:length(dirN)
    cd(dirN{iSess})
    basepath = cd; basename = bz_BasenameFromBasepath(cd);
    
    if exist([basename '.waveforms.mat'],'file')
        load([basename '.waveforms.mat'])
    else
        [waveforms] = pullWideSpike(basepath);
    end
    
    [~,~,aacs] = splitCellTypes(basepath);
    
    [spkwidth] = getSpkWidth(basepath, waveforms,'units',aacs);

    allSpkWidth = [allSpkWidth , spkwidth];
    allSess     = [allSess, repmat(iSess,1,length(aacs))];
    allAACIndx  = [allAACIndx, aacs];
end