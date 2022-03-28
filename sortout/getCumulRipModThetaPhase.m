
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
    ph_mod_5_10Hz = load([basename '.ph_mod_5_10Hz.mat']);
    ph_mod_39_50Hz = load([basename '.ph_mod_39_50Hz.mat']);
    load([basename '.ripmod.mat'])
    load([basename '.cell_metrics.cellinfo.mat'])
    load([basename '_celltypes.mat']);
    
    if ~isempty(aacs)
        
        cumul_gammamod_aac      = [cumul_gammamod_aac, ph_mod_39_50Hz.ph_mod.ph_pref(aacs)];
        cumul_thetaphase_aac    = [cumul_thetaphase_aac, ph_mod_5_10Hz.ph_mod.ph_pref(aacs)];
        cumul_ripmod_aac        = [cumul_ripmod_aac, ripmod.mod(aacs)];
        
        for iAAC = aacs
            IDCt = IDCt +1;
            cumul_ID            = [cumul_ID {[num2str(iSess) '_' num2str(iAAC)]}];
        end
        
        cumul_gammamod_pyr      = [cumul_gammamod_pyr, ph_mod_39_50Hz.ph_mod.ph_pref(pyrs)];
        cumul_thetaphase_pyr    = [cumul_thetaphase_pyr, ph_mod_5_10Hz.ph_mod.ph_pref(pyrs)];
        cumul_ripmod_pyr        = [cumul_ripmod_pyr, ripmod.mod(pyrs)];
        
        cumul_gammamod_ints     = [cumul_gammamod_ints, ph_mod_39_50Hz.ph_mod.ph_pref(ints)];
        cumul_thetaphase_ints   = [cumul_thetaphase_ints, ph_mod_5_10Hz.ph_mod.ph_pref(ints)];
        cumul_ripmod_ints       = [cumul_ripmod_ints, ripmod.mod(ints)];
    end
end
