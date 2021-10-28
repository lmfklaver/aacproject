
for iSess = 17
    
    cd(dirN{iSess})
    basepath = cd; 
    basename = bz_BasenameFromBasepath(cd);
    
    clear optoStim
    load([basename '.optostim.manipulation.mat'])
    
[zeta] = runZeta(basepath,optoStim.timestamps(:,1),'saveMat',true)
end
