launchDirNforAACSessions

%%
for iSess =[3:5, 8, 9, 15:length(dirN)]%  [1:5, 8:length(dirN)]
    
    cd(dirN{iSess})
    
basepath = cd; basename = bz_BasenameFromBasepath(cd);
ConvertPhyKS2toCluResSpk(basepath,basename,[])

end
