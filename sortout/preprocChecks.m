% % check  ait en optostim

for iSess = 1:length(dirN)
    cd(dirN{iSess})
    
basepath = cd;, basename = bz_BasenameFromBasepath(basepath)
ait = LoadEvents([basename '.evt.ait']);
load([basename '.optoStim.manipulation.mat'])
length(ait.time)/2 == length(optoStim.timestamps)

end
