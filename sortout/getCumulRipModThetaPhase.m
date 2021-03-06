
% % % % % % % % % % % % % % % % % % % % % % % %
% % % Get cumulative values AACs for
% % % theta phase and ripple modulation
% % % % % % % % % % % % % % % % % % % % % % % %

plotCount = 0;
figure

cumul_thetaphase_aac = [];
cumul_gammamod_aac = [];
cumul_ripmod_aac = [];
cumul_thetaphase_pyr = [];
cumul_gammamod_pyr = [];
cumul_ripmod_pyr = [];
cumul_thetaphase_ints = [];
cumul_gammamod_ints = [];
cumul_ripmod_ints = [];

cumul_ID=[];
IDCt = 0;

for iSess = sessions
    clear optmod ph_mod STP cell_metrics
    
    cd(dirN{iSess})
    basepath    = cd;
    basename    = bz_BasenameFromBasepath(basepath);
    
    
    load([basename '.STP.mat'])
    load([basename '.ph_mod.mat'])
    load([basename '.cell_metrics.cellinfo.mat'])
    
    [pyrs, ints, aacs] = splitCellTypes(basepath);
    
    if ~isempty(aacs)
        
        cumul_gammamod_aac      = [cumul_gammamod_aac, ph_mod.ph_pref_gam(aacs)];
        cumul_thetaphase_aac    = [cumul_thetaphase_aac, ph_mod.ph_pref_theta(aacs)];
        cumul_ripmod_aac        = [cumul_ripmod_aac, ph_mod.ripmod.mod(aacs)'];
        
        for iAAC = aacs
            IDCt = IDCt +1;
            cumul_ID            = [cumul_ID {[num2str(iSess) '_' num2str(iAAC)]}];
        end
        
        cumul_gammamod_pyr      = [cumul_gammamod_pyr, ph_mod.ph_pref_gam(pyrs)];
        cumul_thetaphase_pyr    = [cumul_thetaphase_pyr, ph_mod.ph_pref_theta(pyrs)];
        cumul_ripmod_pyr        = [cumul_ripmod_pyr, ph_mod.ripmod.mod(pyrs)'];
        
        cumul_gammamod_ints     = [cumul_gammamod_ints, ph_mod.ph_pref_gam(ints)];
        cumul_thetaphase_ints   = [cumul_thetaphase_ints, ph_mod.ph_pref_theta(ints)];
        cumul_ripmod_ints       = [cumul_ripmod_ints, ph_mod.ripmod.mod(ints)'];
    end
end
