function [cumul_ripmod] = getCumulRipMod(basepath, epochs, sessions, varargin)

figure,

cumul_ripmod_aac = [];
cumul_ripmod_pyr = [];
cumul_ripmod_ints = [];
cumul_ID=[];
cumulclu_ID=[];
IDCt = 0;

for iSess = sessions
    cd(dirN{iSess})
    basepath    = cd;
    basename    = bz_BasenameFromBasepath(basepath);
    load([basename '.ripmod.mat'])
    load([basename '.ripple_ccg.mat'])
    load([basename '.spikes.cellinfo.mat'])
    [pyrs, ints, aacs] = splitCellTypes(basepath);
     if ~isempty(aacs)
          cumul_ripmod.cumul_ripmod_aac        = [cumul_ripmod_aac, ph_mod.ripmod.mod(aacs)'];
                for iAAC = aacs
                IDCt = IDCt +1;
                cumul_ID                = [cumul_ID {[num2str(iSess) '_' num2str(iAAC)]}];
                cumulclu_ID             = [cumulclu_ID {[num2str(iSess) '_' num2str(spikes.cluID(iAAC))]}]
                end
          cumul_ripmod.cumul_ripmod_pyr        = [cumul_ripmod_pyr, ph_mod.ripmod.mod(pyrs)'];
          cumul_ripmod.cumul_ripmod_ints       = [cumul_ripmod_ints, ph_mod.ripmod.mod(ints)'];
     end
 if saveMat
save([basename '.cumul_ripmod.mat'],'cumul_ripmod')
end
end
