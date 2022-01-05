
% % % % % % % % % % % % % % % % % % % % % % % %
% % % Get cumulative values AACs for
% % % theta phase and ripple modulation
% % % % % % % % % % % % % % % % % % % % % % % %

% NB: requires a sessions variable to work!

plotCount = 0;
figure

cumul_thetaphase_aac = []; cumul_gammamod_aac = []; cumul_ripmod_aac = [];
cumul_thetaphase_pyr = []; cumul_gammamod_pyr = []; cumul_ripmod_pyr = [];
cumul_thetaphase_ints = []; cumul_gammamod_ints = []; cumul_ripmod_ints = [];

cumul_ID=[];
IDCt = 0;

pref_gamma = ph_mod.ph_pref_gam;
pref_theta = ph_mod.ph_pref_theta;
ripmod = ph_mod.ripmod.mod;

for iSess = sessions
    clear optmod ph_mod STP cell_metrics
    
    cd(dirN{iSess})
    basepath    = cd;
    basename    = bz_BasenameFromBasepath(basepath);
    
    
    load([basename '.STP.mat']) % made in getShort
    load([basename '.ph_mod.mat'])
    load([basename '.cell_metrics.cellinfo.mat'])
   
    
    [pyrs, ints, aacs] = splitCellTypes(basepath);
    
   
    if ~isempty(aacs)
        
        cumul_gammamod_aac      = [cumul_gammamod_aac, pref_gamma(aacs)];
        cumul_thetaphase_aac    = [cumul_thetaphase_aac, pref_theta(aacs)];
        cumul_ripmod_aac        = [cumul_ripmod_aac, ripmod(aacs)'];
        
        for iAAC = aacs
            IDCt = IDCt +1;
            cumul_ID            = [cumul_ID {[num2str(iSess) '_' num2str(iAAC)]}];
        end
        
        cumul_gammamod_pyr      = [cumul_gammamod_pyr, pref_gamma(pyrs)];
        cumul_thetaphase_pyr    = [cumul_thetaphase_pyr, pref_theta(pyrs)];
        cumul_ripmod_pyr        = [cumul_ripmod_pyr, ripmod(pyrs)'];
        
        cumul_gammamod_ints     = [cumul_gammamod_ints, pref_gamma(ints)];
        cumul_thetaphase_ints   = [cumul_thetaphase_ints, pref_theta(ints)];
        cumul_ripmod_ints       = [cumul_ripmod_ints, ripmod(ints)'];
    end
end
