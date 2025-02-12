
% % % % % % % % % % % % % % % % % % % % % % % %
% % % Get cumulative values AACs for
% % % theta phase and ripple modulation
% % % % % % % % % % % % % % % % % % % % % % % %

plotCount = 0;

cumul_thetaphase_aac = [];
cumul_ripphase_aac = [];
cumul_gammamod_aac = [];
cumul_ripmod_aac = [];
cumul_ripmodCE_aac      = [];
cumul_ripmodMAX_aac     = [];
cumul_ripmodMIN_aac     = [];
cumul_ripmodBL_aac      = [];
cumul_ripmodPETH_aac    = [];
cumul_thetaphase_pyr = [];
cumul_ripphase_pyr = [];
cumul_gammamod_pyr = [];
cumul_ripmod_pyr = [];
cumul_ripmodCE_pyr      = [];
cumul_ripmodMAX_pyr     = [];
cumul_ripmodMIN_pyr     = [];
cumul_ripmodBL_pyr      = [];
cumul_ripmodPETH_pyr    = [];
cumul_thetaphase_ints = [];
cumul_ripphase_ints = [];
cumul_gammamod_ints = [];
cumul_ripmod_ints = [];
cumul_ripmodCE_ints      = [];
cumul_ripmodMAX_ints     = [];
cumul_ripmodMIN_ints     = [];
cumul_ripmodBL_ints      = [];
cumul_ripmodPETH_ints    = [];
cumul_thetamod_aac =[];
cumul_thetamod_pyr=[];
cumul_thetamod_ints=[];
cumul_ID=[];
cumul_numpremono_aac=[];
cumul_numpremono_aacnorm=[];
cumul_ripmodpremono_aac=[];
cumul_ripmodpremonoCE_aac   = [];
cumul_ripmodpremonoCECONTROL_aac   = [];
cumul_ripmodpremonoMAX_aac  = [];
cumul_ripmodpremonoMIN_aac  = [];
cumul_ripmodpremonoBL_aac   = [];
cumul_ripmodpremonoPETH_aac = [];
cumul_FRpremono=[];
cumul_burstpremono=[];
cumul_thetapremono=[];
cumul_ripzeta_aac=[];
cumul_AACripplepeth.rate=[];
cumul_AACripplepeth5s.rate=[];
cumul_AACthetapeth.rate=[];
cumul_AACThetaHist=[];
cumul_meanpremonotransprob=[];
cumul_maxpremonotransprob=[];
cumul_ripmodpremonoMODIFIED=[];
cumul_ripratepowcorrpremono = [];
cumul_ripspikephasemagpremono   = [];
cumul_ripspikephaseanglepremono = [];
cumul_ratepowercorr_aac = [];
cumul_spikephaseangle_aac =  [];
cumul_spikephasemag_aac =  [];
cumul_ratepowercorr_pyr =  [];
cumul_spikephaseangle_pyr =  [];
cumul_spikephasemag_pyr =  [];
cumul_FR_aac=[];
cumul_RipTTFS_aac=[];
cumul_RipSAP_aac=[];
MasAACPairSpkC=[];
cumul_RipCCG_aac=[];
cumul_ripmodDpremonoCE_aac   = [];
cumul_ripmodSpremonoCE_aac   = [];
cumul_waveform_aac=[];
cumul_pos_aac=[];
cumul_burstiness_aac=[];
cumul_isi_aac=[];
cumul_ripisi_aac=[];
cumul_abratio_aac=[];
cumul_spikeamp_aac=[];
cumul_phy_amp_aac=[];
cumul_cv2_aac=[];
cumul_t2p_aac=[];
cumul_ripspkmean_aac    = [];
cumul_ripspkmax_aac     = [];
cumul_acg_aac=[];
cumul_thetap_aac=[];
cumul_thetar_aac=[];
cumul_thetam_aac=[];
cumul_ripplep_aac=[];
cumul_rippler_aac=[];
cumul_ripplem_aac=[];
cumul_ripccg_aac=[];
cumul_thetaanglepremono=[];
cumul_thetaanglepremononorm=[];
cumul_thetamagpremono=[];
cumul_ripmodCESig_aac=[];
cumul_ripmodNoPre_aac=[];
cumul_ripplep_premono=[];
cumul_rippler_premono=[];
cumul_ripplem_premono=[];
cumul_AACSpk=[];
cumul_AACPreSynSpk=[];
cumul_AACNonConSpk=[];
for iSess = sessions
    clear optmod ph_mod STP cell_metrics
    numpremono=[];
    numpremononorm=[];
    ripmodpremonoCE=[];
    ripmodpremonoCEcontrol=[];
    ripmodDpremonoCE=[];
    ripmodSpremonoCE=[];
    ripmodpremono=[];
    ripmodpremonoMAX=[];
    ripmodpremonoMIN=[];
    ripmodpremonoBL=[];
    ripmodpremonoPETH=[];
    ripmodpremonoMODIFIED=[];
    ripmodNoPre=[];
    ripratepowcorrpremono=[];
    ripspikephasemagpremono=[];
    ripspikephaseanglepremono=[];
    FRpremono=[];
    burstpremono=[];
    thetapremono=[];
    AACripplepeth=[];
    AACripplepeth5s=[];
    AACthetapeth=[];
    AACThetaHist=[];
    meanpremonotransprob=[];
    maxpremonotransprob=[];
    thetaanglepremono=[];
    thetaanglepremononorm=[];
    thetamagpremono=[];
    ripplep_premono=[];
    rippler_premono=[];
    ripplem_premono=[];
    AACPairSpkc=[];
    AACSpk=[];
    AACPreSynSpk=[];
    AACNonConSpk=[];
    cd(dirN{iSess})
    basepath    = cd;
    basename    = bz_BasenameFromBasepath(basepath);
    if isfile([basename '.ripspikes.allripinstim.analysis.mat'])
            load([basename '.ripzeta.stats.mat']);
    load([basename '.ripspikes.allripinstim.analysis.mat']);
    load([basename '.ripspiketime.analysis.mat'])
    else 
    end
    load([basename '.SpikeLFPCouplingGdEps1F.mat'])
    ph_mod_5_10Hz = load([basename '.ph_mod_5_10Hz.mat']);
    ph_mod_39_50Hz = load([basename '.ph_mod_39_50Hz.mat']);
    ph_mod_120_250Hz = load([basename '.ph_mod_120_250Hz.mat']);
    load([basename '.ripmod.mat']);
    %load([basename '.cell_metrics.cellinfo_RipManip.mat']); %Uncomment and
    %comment next line to build NoRipCellExplorer for control values
    load([basename '.cell_metrics.cellinfo.mat']);
    load([basename '_celltypes.mat']);
    load([basename '.thetamodulationindex.mat']);
    load([basename '.burstMizuseki.analysis.mat']);
    load([basename '.ripzeta.stats.mat']);
    load([basename '.ripple_ccg.mat'])
    load([basename '.spikes.cellinfo.mat'])
    load([basename '.ripspikes.NoStim.analysis.mat'])
    load([basename '.ripSTA.firstripspikeonly.analysis.mat'])
    rippeth=load([basename '.ripplepeth.analysis.mat']);
    rippeth5s=load([basename '.ripplepeth5s.analysis.mat']);
    thetapeth=load([basename '.thetapeth.analysis.mat']);
    load([basename '.theta_6-12.PhaseLockingData.cellinfo.mat']);
    load([basename '.ripple_100-250.PhaseLockingData.cellinfo.mat']);
    if ~isempty(aacs)
        thetapref=ph_mod_5_10Hz.ph_mod.ph_pref-pi;
        switchind=find(thetapref<-pi);
        thetapref(switchind)=thetapref(switchind)+2*pi;
        switchind=[];
        rippref=ph_mod_120_250Hz.ph_mod.ph_pref-pi;
        switchind=find(rippref<-pi);
        rippref(switchind)=rippref(switchind)+2*pi;
        cumul_acg_aac      = [cumul_acg_aac, cell_metrics.acg.narrow(:,aacs)];
        cumul_gammamod_aac      = [cumul_gammamod_aac, ph_mod_39_50Hz.ph_mod.ph_pref(aacs)];
        cumul_ripphase_aac      = [cumul_ripphase_aac, rippref(aacs)];
        %cumul_thetaphase_aac    = [cumul_thetaphase_aac, thetapref(aacs)];
        cumul_ripmodCE_aac      = [cumul_ripmodCE_aac,  cell_metrics.ripples_modulationIndex(aacs)];
        cumul_ripmodCESig_aac   = [cumul_ripmodCESig_aac cell_metrics.ripples_modulationSignificanceLevel(aacs)];
        cumul_ripmod_aac        = [cumul_ripmod_aac,  ripmod.mod(aacs)];
        cumul_ripmodMAX_aac     = [cumul_ripmodMAX_aac, ripmod.max(aacs)];
        cumul_ripmodMIN_aac     = [cumul_ripmodMIN_aac, ripmod.min(aacs)];
        cumul_ripmodBL_aac      = [cumul_ripmodBL_aac, ripmod.meanbaseline(aacs)];
        cumul_ripmodPETH_aac    = [cumul_ripmodPETH_aac, ripmod.modPETH(aacs)];
        cumul_ripzeta_aac       = [cumul_ripzeta_aac, zeta.P(aacs)];
        cumul_ratepowercorr_aac = [cumul_ratepowercorr_aac, SpikeLFPCouplingGdEps.cell.ratepowercorr(aacs,:,1)'];
        cumul_spikephaseangle_aac = [cumul_spikephaseangle_aac ,SpikeLFPCouplingGdEps.cell.spikephaseangle(aacs,:,1)'];
        cumul_spikephasemag_aac = [cumul_spikephasemag_aac ,SpikeLFPCouplingGdEps.cell.spikephasemag(aacs,:,1)'];
        cumul_thetamod_aac      = [cumul_thetamod_aac, cell_metrics.thetaModulationIndex(aacs)];
        cumul_FR_aac            = [cumul_FR_aac, cell_metrics.firingRate(aacs)];
        cumul_RipTTFS_aac       = [cumul_RipTTFS_aac, ripspiketime.RipSpkTimeStrtMean(aacs)];
        cumul_RipSAP_aac        = [cumul_RipSAP_aac, ripspiketime.RipSpkTimePeakMean(aacs)];
        cumul_RipCCG_aac        = [cumul_RipCCG_aac, squeeze(ripple_ccg.ccg(:,end,aacs))];
        cumul_waveform_aac      = [cumul_waveform_aac; cell2mat(cell_metrics.waveforms.filt(:, aacs)')];
        cumul_pos_aac           = [cumul_pos_aac, cell_metrics.deepSuperficialDistance(aacs)];
        cumul_burstiness_aac    = [cumul_burstiness_aac, cell_metrics.burstIndex_Mizuseki2012(aacs)];
        cumul_isi_aac           = [cumul_isi_aac, cell_metrics.firingRateISI(aacs)];
        cumul_abratio_aac       = [cumul_abratio_aac, cell_metrics.ab_ratio(aacs)];
        cumul_spikeamp_aac      = [cumul_spikeamp_aac, cell_metrics.peakVoltage(aacs)];
        cumul_phy_amp_aac       = [cumul_phy_amp_aac, cell_metrics.phy_amp(aacs)];
        cumul_cv2_aac           = [cumul_cv2_aac, cell_metrics.cv2(aacs)];
        cumul_t2p_aac           = [cumul_t2p_aac, cell_metrics.troughToPeak(aacs)];
        cumul_ripisi_aac        = [cumul_ripisi_aac, ripSTA.ripISIMean(aacs)];
        ripspikesNoStim.numSpkperRip_OFF(ripspikesNoStim.numSpkperRip_OFF==0)=nan;
        cumul_ripspkmean_aac    = [cumul_ripspkmean_aac nanmean(ripspikesNoStim.numSpkperRip_OFF(aacs,:),2)'];
        cumul_ripspkmax_aac     = [cumul_ripspkmax_aac max(ripspikesNoStim.numSpkperRip_OFF(aacs,:),[],2)'];
        cumul_thetap_aac        = [cumul_thetap_aac thetaMod.phasestats.p(aacs)];
        cumul_thetar_aac        = [cumul_thetar_aac thetaMod.phasestats.r(aacs)];
        cumul_thetam_aac        = [cumul_thetam_aac thetaMod.phasestats.m(aacs)];
        cumul_ripplep_aac        = [cumul_ripplep_aac rippleMod.phasestats.p(aacs)];
        cumul_rippler_aac        = [cumul_rippler_aac rippleMod.phasestats.r(aacs)];
        cumul_ripplem_aac        = [cumul_ripplem_aac rippleMod.phasestats.m(aacs)];
        for iAAC = aacs;
            cumul_ID            = [cumul_ID {[num2str(iSess) '_' num2str(iAAC)]}];
        end
        
        cumul_gammamod_pyr      = [cumul_gammamod_pyr, ph_mod_39_50Hz.ph_mod.ph_pref(pyrs)];
        cumul_thetaphase_pyr    = [cumul_thetaphase_pyr, thetapref(pyrs)];
        cumul_ripphase_pyr      = [cumul_ripphase_pyr, rippref(pyrs)];        
        cumul_ripmodCE_pyr      = [cumul_ripmodCE_pyr,  cell_metrics.ripples_modulationIndex(pyrs)];
        cumul_ripmod_pyr        = [cumul_ripmod_pyr,  ripmod.mod(pyrs)];
        cumul_ripmodMAX_pyr     = [cumul_ripmodMAX_pyr, ripmod.max(pyrs)];
        cumul_ripmodMIN_pyr     = [cumul_ripmodMIN_pyr, ripmod.min(pyrs)];
        cumul_ripmodBL_pyr      = [cumul_ripmodBL_pyr, ripmod.meanbaseline(pyrs)];
        cumul_ripmodPETH_pyr    = [cumul_ripmodPETH_pyr, ripmod.modPETH(pyrs)];
        cumul_thetamod_pyr      = [cumul_thetamod_pyr, thetaModulationIndex(pyrs)];
        cumul_ratepowercorr_pyr = [cumul_ratepowercorr_pyr, SpikeLFPCouplingGdEps.cell.ratepowercorr(pyrs,:,1)'];
        cumul_spikephaseangle_pyr = [cumul_spikephaseangle_pyr ,SpikeLFPCouplingGdEps.cell.spikephaseangle(pyrs,:,1)'];
        cumul_spikephasemag_pyr = [cumul_spikephasemag_pyr ,SpikeLFPCouplingGdEps.cell.spikephasemag(pyrs,:,1)'];
        
        cumul_gammamod_ints     = [cumul_gammamod_ints, ph_mod_39_50Hz.ph_mod.ph_pref(ints)];
        cumul_thetaphase_ints   = [cumul_thetaphase_ints, thetapref(ints)];
        cumul_ripphase_ints      = [cumul_ripphase_ints, rippref(ints)];        

        cumul_ripmodCE_ints      = [cumul_ripmodCE_ints,  cell_metrics.ripples_modulationIndex(ints)];
        cumul_ripmod_ints        = [cumul_ripmod_ints,  ripmod.mod(ints)];
        cumul_ripmodMAX_ints     = [cumul_ripmodMAX_ints, ripmod.max(ints)];
        cumul_ripmodMIN_ints     = [cumul_ripmodMIN_ints, ripmod.min(ints)];
        cumul_ripmodBL_ints      = [cumul_ripmodBL_ints, ripmod.meanbaseline(ints)];
        cumul_ripmodPETH_ints    = [cumul_ripmodPETH_ints, ripmod.modPETH(ints)];
        cumul_thetamod_ints      = [cumul_thetamod_ints, thetaModulationIndex(ints)];
        for nAAC = 1:length(aacs);
            transprobs=cell_metrics.putativeConnections.excitatoryTransProb(find(cell_metrics.putativeConnections.excitatory(:,2)==aacs(nAAC)));
            meanpremonotransprob(nAAC)=nanmean(transprobs);
            if isempty(transprobs);
                maxpremonotransprob(nAAC)=NaN;
            else
            maxpremonotransprob(nAAC)=max(transprobs);
            end
            presynmono=cell_metrics.putativeConnections.excitatory((find(cell_metrics.putativeConnections.excitatory(:,2) == aacs(nAAC))));
            presynmono=presynmono(ismember(presynmono,pyrs));
            numpremono(nAAC)=length(presynmono);
            numpremononorm(nAAC)=length(presynmono)/length(pyrs);
            Dindices=find(contains(cell_metrics.deepSuperficial,'Deep'));
            Sindices=find(contains(cell_metrics.deepSuperficial,'Superficial'));
            Dpresynmono=presynmono(ismember(presynmono,Dindices));
            Spresynmono=presynmono(ismember(presynmono,Sindices));
            ripmodpremonoCE(nAAC)=nanmean(cell_metrics.ripples_modulationIndex(presynmono));
            if ~isempty(presynmono)
                ripmodNoPre(nAAC)=nanmean(cell_metrics.ripples_modulationIndex(pyrs(~ismember(pyrs,presynmono))));
                else
                ripmodNoPre(nAAC)=nanmean(cell_metrics.ripples_modulationIndex(pyrs));
            end
            % Step 1: Identify pyramidal cells not in presynmono
            non_presynmono_pyrs = pyrs(~ismember(pyrs, presynmono));
            % Step 2: Randomly select the same number of cells as in presynmono
            num_presynmono = length(presynmono);
            random_indices = randperm(length(non_presynmono_pyrs), num_presynmono);
            random_subset = non_presynmono_pyrs(random_indices);
            ripmodpremonoCEcontrol(nAAC) = nanmean(cell_metrics.ripples_modulationIndex(random_subset));
            ripmodDpremonoCE(nAAC)=nanmean(cell_metrics.ripples_modulationIndex(Dpresynmono));
            ripmodSpremonoCE(nAAC)=nanmean(cell_metrics.ripples_modulationIndex(Spresynmono));
            ripmodpremono(nAAC)=nanmean(ripmod.mod(presynmono(ismember(presynmono,pyrs))));
            ripmodpremonoMAX(nAAC)=nanmean(ripmod.max(presynmono(ismember(presynmono,pyrs))));
            ripmodpremonoMIN(nAAC)=nanmean(ripmod.min(presynmono(ismember(presynmono,pyrs))));
            ripmodpremonoBL(nAAC)=nanmean(ripmod.meanbaseline(presynmono(ismember(presynmono,pyrs))));
            ripmodpremonoPETH(nAAC)=nanmean(ripmod.modPETH(presynmono(ismember(presynmono,pyrs))));
            ripmodpremonoMODIFIED(nAAC)=nanmean((ripmod.max(presynmono(ismember(presynmono,pyrs)))-ripmod.min(presynmono(ismember(presynmono,pyrs))))/ripmod.meanbaseline(presynmono(ismember(presynmono,pyrs))));
            ripratepowcorrpremono(nAAC)=nanmean(SpikeLFPCouplingGdEps.cell.ratepowercorr(presynmono(ismember(presynmono,pyrs)),:,1));
            ripspikephasemagpremono(nAAC)=nanmean(SpikeLFPCouplingGdEps.cell.spikephasemag(presynmono(ismember(presynmono,pyrs)),:,1));
            ripspikephaseanglepremono(nAAC)=nanmean(SpikeLFPCouplingGdEps.cell.spikephaseangle(presynmono(ismember(presynmono,pyrs)),:,1));
            ripplep_premono(nAAC)    =nanmean(rippleMod.phasestats.p(presynmono(ismember(presynmono,pyrs))));
            rippler_premono(nAAC)    =nanmean(rippleMod.phasestats.r(presynmono(ismember(presynmono,pyrs))));
            ripplem_premono(nAAC)    =nanmean(rippleMod.phasestats.m(presynmono(ismember(presynmono,pyrs))));
            FRpremono(nAAC)=mean(cell_metrics.firingRate(presynmono(ismember(presynmono,pyrs))));
            burstpremono(nAAC)=mean(burstIndex.in(presynmono(ismember(presynmono,pyrs))));
            thetapremono(nAAC)=mean(thetaModulationIndex(presynmono(ismember(presynmono,pyrs))));
            pyrpresynmono=presynmono(ismember(presynmono,pyrs));
            sigthetapresynmono=pyrpresynmono(find(thetaMod.phasestats.p(pyrpresynmono)<0.05));
            thetaanglepremono(nAAC)=mean(thetaMod.phasestats.m(sigthetapresynmono));
            thetaanglepremononorm(nAAC)=mean(thetaMod.phasestats.m(sigthetapresynmono).*thetaMod.phasestats.r(sigthetapresynmono));
            thetamagpremono(nAAC)=mean(thetaMod.phasestats.r(sigthetapresynmono));
            AACripplepeth(nAAC,:)=rippeth.ripplepeth.rate(aacs(nAAC),:);
            AACripplepeth5s(nAAC,:)=rippeth5s.ripplepeth5s.rate(aacs(nAAC),:);
            AACthetapeth(nAAC,:)=thetapeth.thetapeth.rate(aacs(nAAC),:);
            AACPairSpk=[];
            AACPairSpk(1,:)=ripspikes.numSpkperRip_OFF(aacs(nAAC),:);
            AACPairSpkc{nAAC}=[AACPairSpk; ripspikes.numSpkperRip_OFF(presynmono,:)];
            AACSpk{nAAC}=spikes.times(aacs(nAAC));
            AACPreSynSpk{nAAC}=spikes.times(:,presynmono);
            AACNonConSpk{nAAC}=spikes.times(:,non_presynmono_pyrs);
        end
        ripmodpremono(isnan(ripmodpremono))=0;
        ripmodpremonoCE(isnan(ripmodpremonoCE))=0;
        ripmodpremonoCEcontrol(isnan(ripmodpremonoCEcontrol))=0;
        ripmodDpremonoCE(isnan(ripmodDpremonoCE))=0;
        ripmodSpremonoCE(isnan(ripmodSpremonoCE))=0;
        ripmodpremonoMAX(isnan(ripmodpremonoMAX))=0;
        ripmodpremonoMIN(isnan(ripmodpremonoMIN))=0;
        ripmodpremonoBL(isnan(ripmodpremonoBL))=0;
        ripmodpremonoPETH(isnan(ripmodpremonoPETH))=0;
        ripratepowcorrpremono(isnan(ripratepowcorrpremono))=0;
        ripspikephasemagpremono(isnan(ripspikephasemagpremono))=0;
        ripspikephaseanglepremono(isnan(ripspikephaseanglepremono))=0;
        cumul_numpremono_aac        = [cumul_numpremono_aac, numpremono];
        cumul_numpremono_aacnorm    = [cumul_numpremono_aacnorm, numpremononorm];
        cumul_ripmodpremono_aac     = [cumul_ripmodpremono_aac, ripmodpremono];
        cumul_ripmodpremonoCE_aac   = [cumul_ripmodpremonoCE_aac, ripmodpremonoCE];
        cumul_ripmodpremonoCECONTROL_aac = [cumul_ripmodpremonoCECONTROL_aac, ripmodpremonoCEcontrol];
        cumul_ripmodNoPre_aac       = [cumul_ripmodNoPre_aac, ripmodNoPre];
        cumul_ripmodDpremonoCE_aac   = [cumul_ripmodDpremonoCE_aac, ripmodDpremonoCE];
        cumul_ripmodSpremonoCE_aac   = [cumul_ripmodSpremonoCE_aac, ripmodSpremonoCE];
        cumul_ripmodpremonoMAX_aac  = [cumul_ripmodpremonoMAX_aac, ripmodpremonoMAX];
        cumul_ripmodpremonoMIN_aac  = [cumul_ripmodpremonoMIN_aac, ripmodpremonoMIN];
        cumul_ripmodpremonoBL_aac   = [cumul_ripmodpremonoBL_aac, ripmodpremonoBL];
        cumul_ripmodpremonoPETH_aac = [cumul_ripmodpremonoPETH_aac, ripmodpremonoPETH];
        cumul_ripmodpremonoMODIFIED = [cumul_ripmodpremonoMODIFIED, ripmodpremonoMODIFIED];
        cumul_ripplep_premono    = [cumul_ripplep_premono ripplep_premono];
        cumul_rippler_premono    = [cumul_rippler_premono rippler_premono];
        cumul_ripplem_premono    = [cumul_ripplem_premono ripplem_premono];
        cumul_ripratepowcorrpremono = [cumul_ripratepowcorrpremono, ripratepowcorrpremono];
        cumul_ripspikephasemagpremono   = [cumul_ripspikephasemagpremono, ripspikephasemagpremono];
        cumul_ripspikephaseanglepremono = [cumul_ripspikephaseanglepremono, ripspikephaseanglepremono];   
        cumul_FRpremono         = [cumul_FRpremono, FRpremono];
        cumul_burstpremono      = [cumul_burstpremono, burstpremono];
        cumul_thetapremono      = [cumul_thetapremono, thetapremono];
        cumul_thetaanglepremono=[cumul_thetaanglepremono thetaanglepremono];
        cumul_thetaanglepremononorm=[cumul_thetaanglepremononorm thetaanglepremononorm];
        cumul_thetamagpremono=[cumul_thetamagpremono thetamagpremono];
        cumul_AACripplepeth.rate     = [cumul_AACripplepeth.rate; AACripplepeth];
        cumul_AACripplepeth5s.rate     = [cumul_AACripplepeth5s.rate; AACripplepeth5s];
        cumul_AACthetapeth.rate     = [cumul_AACthetapeth.rate; AACthetapeth];
        cumul_AACripplepeth.timeEdges= rippeth.ripplepeth.timeEdges;
        cumul_AACripplepeth5s.timeEdges= rippeth5s.ripplepeth5s.timeEdges;
        cumul_AACThetaHist      =[cumul_AACThetaHist; AACThetaHist];
        cumul_meanpremonotransprob  =[cumul_meanpremonotransprob, meanpremonotransprob];
        cumul_maxpremonotransprob  =[cumul_maxpremonotransprob, maxpremonotransprob];
        
    MasAACPairSpkC=[MasAACPairSpkC AACPairSpkc];
            cumul_AACSpk=[cumul_AACSpk AACSpk];
            cumul_AACPreSynSpk=[cumul_AACPreSynSpk AACPreSynSpk];
            cumul_AACNonConSpk=[cumul_AACNonConSpk AACNonConSpk];
    end

end
