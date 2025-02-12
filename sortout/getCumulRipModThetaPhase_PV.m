
% % % % % % % % % % % % % % % % % % % % % % % %
% % % Get cumulative values pvs for
% % % theta phase and ripple modulation
% % % % % % % % % % % % % % % % % % % % % % % %
sessions=[24:34]
plotCount = 0;

cumul_thetaphase_pv = [];
cumul_ripphase_pv = [];
cumul_gammamod_pv = [];
cumul_thetamod_pv = [];
cumul_ripmod_pv = [];
cumul_ripmodCE_pv      = [];
cumul_ripmodMAX_pv     = [];
cumul_ripmodMIN_pv     = [];
cumul_ripmodBL_pv      = [];
cumul_ripmodPETH_pv    = [];
cumul_thetaphase_pyr_pvses = [];
cumul_ripphase_pyr_pvses = [];
cumul_gammamod_pyr_pvses = [];
cumul_ripmod_pyr_pvses = [];
cumul_ripmodCE_pyr_pvses      = [];
cumul_ripmodMAX_pyr_pvses     = [];
cumul_ripmodMIN_pyr_pvses     = [];
cumul_ripmodBL_pyr_pvses      = [];
cumul_ripmodPETH_pyr_pvses    = [];
cumul_ID=[];
cumul_numpremono_pv=[];
cumul_numpremono_pvnorm=[];
cumul_ripmodpremono_pv=[];
cumul_ripmodpremonoCE_pv   = [];
cumul_ripmodpremonoMAX_pv  = [];
cumul_ripmodpremonoMIN_pv  = [];
cumul_ripmodpremonoBL_pv   = [];
cumul_ripmodpremonoPETH_pv = [];
cumul_FRpremono=[];
cumul_burstpremono=[];
cumul_thetapremono=[];
cumul_ripzeta_pv=[];
cumul_pvripplepeth.rate=[];
cumul_pvThetaHist=[];
cumul_meanpremonotransprob=[];
cumul_maxpremonotransprob=[];
cumul_ripmodpremonoMODIFIED=[];
cumul_ripratepowcorrpremono = [];
cumul_ripspikephasemagpremono   = [];
cumul_ripspikephaseanglepremono = [];
cumul_ratepowercorr_pv = [];
cumul_spikephaseangle_pv =  [];
cumul_spikephasemag_pv =  [];
cumul_ratepowercorr_pyr_pvses =  [];
cumul_spikephaseangle_pyr_pvses =  [];
cumul_spikephasemag_pyr_pvses =  [];
cumul_FR_pv=[];
cumul_RipTTFS_pv=[];
cumul_RipSAP_pv=[];
MaspvPairSpkC=[];
cumul_RipCCG_pv=[];
cumul_ripmodDpremonoCE_pv   = [];
cumul_ripmodSpremonoCE_pv   = [];
cumul_waveform_pv=[];
cumul_pos_pv=[];
cumul_burstiness_pv=[];
cumul_isi_pv=[];
cumul_ripisi_pv=[];
cumul_ripspkmean_pv    = [];
cumul_ripspkmax_pv     = [];
cumul_acg_pv=[];
cumul_thetap_pv=[];
cumul_thetar_pv=[];
cumul_thetam_pv=[];
cumul_ripccg_pv=[];
cumul_thetaanglepremono=[];
cumul_thetaanglepremononorm=[];
cumul_thetamagpremono=[];
cumul_ripmodCESig_pv=[];
cumul_ripmodNoPre_pv=[];
cumul_cv2_pv=[];
cumul_t2p_pv=[];
cumul_abratio_pv=[];
cumul_spikeamp_pv=[];
cumul_phy_amp_pv=[];

for iSess = sessions
    clear optmod ph_mod STP cell_metrics
    numpremono=[];
    numpremononorm=[];
    ripmodpremonoCE=[];
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
    pvripplepeth=[];
    pvThetaHist=[];
    meanpremonotransprob=[];
    maxpremonotransprob=[];
    thetaanglepremono=[];
    thetaanglepremononorm=[];
    thetamagpremono=[];
    pvPairSpkc=[];
    cd(dirN{iSess})
    basepath    = cd;
    basename    = bz_BasenameFromBasepath(basepath);
    if isfile([basename '.ripspikes.allripinstim.analysis.mat'])
    load([basename '.ripspikes.allripinstim.analysis.mat']);
    else 
    end
    load([basename '.ripspiketime.analysis.mat'])
    load([basename '.SpikeLFPCouplingGdEps1F.mat'])
    load([basename '.ripmod.mat']);
    load([basename '.cell_metrics.cellinfo.mat']);
    load([basename '_celltypes.mat']);
    load([basename '.ripspikes.NoStim.analysis.mat'])
    load([basename '.ripSTA.firstripspikeonly.analysis.mat'])
    rippeth=load([basename '.ripplepeth.analysis.mat']);
    
    if ~isempty(aacs)
        cumul_acg_pv      = [cumul_acg_pv, cell_metrics.acg.narrow(:,aacs)];
        cumul_ripmodCE_pv      = [cumul_ripmodCE_pv,  cell_metrics.ripples_modulationIndex(aacs)];
        cumul_ripmodCESig_pv   = [cumul_ripmodCESig_pv cell_metrics.ripples_modulationSignificanceLevel(aacs)];
        cumul_ripmod_pv        = [cumul_ripmod_pv,  ripmod.mod(aacs)];
        cumul_ripmodMAX_pv     = [cumul_ripmodMAX_pv, ripmod.max(aacs)];
        cumul_ripmodMIN_pv     = [cumul_ripmodMIN_pv, ripmod.min(aacs)];
        cumul_ripmodBL_pv      = [cumul_ripmodBL_pv, ripmod.meanbaseline(aacs)];
        cumul_ripmodPETH_pv    = [cumul_ripmodPETH_pv, ripmod.modPETH(aacs)];
        cumul_ratepowercorr_pv = [cumul_ratepowercorr_pv, SpikeLFPCouplingGdEps.cell.ratepowercorr(aacs,:,1)'];
        cumul_spikephaseangle_pv = [cumul_spikephaseangle_pv ,SpikeLFPCouplingGdEps.cell.spikephaseangle(aacs,:,1)'];
        cumul_spikephasemag_pv = [cumul_spikephasemag_pv ,SpikeLFPCouplingGdEps.cell.spikephasemag(aacs,:,1)'];
        cumul_thetamod_pv      = [cumul_thetamod_pv, cell_metrics.thetaModulationIndex(aacs)];
        cumul_FR_pv            = [cumul_FR_pv, cell_metrics.firingRate(aacs)];
        cumul_RipTTFS_pv       = [cumul_RipTTFS_pv, ripspiketime.RipSpkTimeStrtMean(aacs)];
        cumul_RipSAP_pv        = [cumul_RipSAP_pv, ripspiketime.RipSpkTimePeakMean(aacs)];
        %cumul_waveform_pv      = [cumul_waveform_pv; cell2mat(cell_metrics.waveforms.filt(:, aacs)')];
        cumul_pos_pv           = [cumul_pos_pv, cell_metrics.deepSuperficialDistance(aacs)];
        cumul_burstiness_pv    = [cumul_burstiness_pv, cell_metrics.burstIndex_Mizuseki2012(aacs)];
        cumul_isi_pv           = [cumul_isi_pv, cell_metrics.firingRateISI(aacs)];
        cumul_ripisi_pv        = [cumul_ripisi_pv, ripSTA.ripISIMean(aacs)];
        ripspikesNoStim.numSpkperRip_OFF(ripspikesNoStim.numSpkperRip_OFF==0)=nan;
        cumul_ripspkmean_pv    = [cumul_ripspkmean_pv nanmean(ripspikesNoStim.numSpkperRip_OFF(aacs,:),2)'];
        cumul_ripspkmax_pv     = [cumul_ripspkmax_pv max(ripspikesNoStim.numSpkperRip_OFF(aacs,:),[],2)'];
        cumul_abratio_pv       = [cumul_abratio_pv, cell_metrics.ab_ratio(aacs)];
        cumul_spikeamp_pv      = [cumul_spikeamp_pv, cell_metrics.peakVoltage(aacs)];
        cumul_cv2_pv           = [cumul_cv2_pv, cell_metrics.cv2(aacs)];
        cumul_t2p_pv           = [cumul_t2p_pv, cell_metrics.troughToPeak(aacs)];
        for ipv = aacs;
            cumul_ID            = [cumul_ID {[num2str(iSess) '_' num2str(ipv)]}];
        end
        
        cumul_ripmodCE_pyr_pvses      = [cumul_ripmodCE_pyr_pvses,  cell_metrics.ripples_modulationIndex(pyrs)];
        cumul_ripmod_pyr_pvses        = [cumul_ripmod_pyr_pvses,  ripmod.mod(pyrs)];
        cumul_ripmodMAX_pyr_pvses     = [cumul_ripmodMAX_pyr_pvses, ripmod.max(pyrs)];
        cumul_ripmodMIN_pyr_pvses     = [cumul_ripmodMIN_pyr_pvses, ripmod.min(pyrs)];
        cumul_ripmodBL_pyr_pvses      = [cumul_ripmodBL_pyr_pvses, ripmod.meanbaseline(pyrs)];
        cumul_ripmodPETH_pyr_pvses    = [cumul_ripmodPETH_pyr_pvses, ripmod.modPETH(pyrs)];
        cumul_ratepowercorr_pyr_pvses = [cumul_ratepowercorr_pyr_pvses, SpikeLFPCouplingGdEps.cell.ratepowercorr(pyrs,:,1)'];
        cumul_spikephaseangle_pyr_pvses = [cumul_spikephaseangle_pyr_pvses ,SpikeLFPCouplingGdEps.cell.spikephaseangle(pyrs,:,1)'];
        cumul_spikephasemag_pyr_pvses = [cumul_spikephasemag_pyr_pvses ,SpikeLFPCouplingGdEps.cell.spikephasemag(pyrs,:,1)'];
        

        for npv = 1:length(aacs);
            transprobs=cell_metrics.putativeConnections.excitatoryTransProb(find(cell_metrics.putativeConnections.excitatory(:,2)==aacs(npv)));
            meanpremonotransprob(npv)=nanmean(transprobs);
            if isempty(transprobs);
                maxpremonotransprob(npv)=NaN;
            else
            maxpremonotransprob(npv)=max(transprobs);
            end
            presynmono=cell_metrics.putativeConnections.excitatory((find(cell_metrics.putativeConnections.excitatory(:,2) == aacs(npv))));
            presynmono=presynmono(ismember(presynmono,pyrs));
            numpremono(npv)=length(presynmono);
            numpremononorm(npv)=length(presynmono)/length(pyrs);
            Dindices=find(contains(cell_metrics.deepSuperficial,'Deep'));
            Sindices=find(contains(cell_metrics.deepSuperficial,'Superficial'));
            Dpresynmono=presynmono(ismember(presynmono,Dindices));
            Spresynmono=presynmono(ismember(presynmono,Sindices));
            ripmodpremonoCE(npv)=nanmean(cell_metrics.ripples_modulationIndex(presynmono));
            if ~isempty(presynmono)
                ripmodNoPre(npv)=nanmean(cell_metrics.ripples_modulationIndex(pyrs(~ismember(pyrs,presynmono))));
                else
                ripmodNoPre(npv)=nanmean(cell_metrics.ripples_modulationIndex(pyrs));
            end
            ripmodDpremonoCE(npv)=nanmean(cell_metrics.ripples_modulationIndex(Dpresynmono));
            ripmodSpremonoCE(npv)=nanmean(cell_metrics.ripples_modulationIndex(Spresynmono));
            ripmodpremono(npv)=nanmean(ripmod.mod(presynmono(ismember(presynmono,pyrs))));
            ripmodpremonoMAX(npv)=nanmean(ripmod.max(presynmono(ismember(presynmono,pyrs))));
            ripmodpremonoMIN(npv)=nanmean(ripmod.min(presynmono(ismember(presynmono,pyrs))));
            ripmodpremonoBL(npv)=nanmean(ripmod.meanbaseline(presynmono(ismember(presynmono,pyrs))));
            ripmodpremonoPETH(npv)=nanmean(ripmod.modPETH(presynmono(ismember(presynmono,pyrs))));
            ripmodpremonoMODIFIED(npv)=nanmean((ripmod.max(presynmono(ismember(presynmono,pyrs)))-ripmod.min(presynmono(ismember(presynmono,pyrs))))/ripmod.meanbaseline(presynmono(ismember(presynmono,pyrs))));
            ripratepowcorrpremono(npv)=nanmean(SpikeLFPCouplingGdEps.cell.ratepowercorr(presynmono(ismember(presynmono,pyrs)),:,1));
            ripspikephasemagpremono(npv)=nanmean(SpikeLFPCouplingGdEps.cell.spikephasemag(presynmono(ismember(presynmono,pyrs)),:,1));
            ripspikephaseanglepremono(npv)=nanmean(SpikeLFPCouplingGdEps.cell.spikephaseangle(presynmono(ismember(presynmono,pyrs)),:,1));
            FRpremono(npv)=mean(cell_metrics.firingRate(presynmono(ismember(presynmono,pyrs))));
            pyrpresynmono=presynmono(ismember(presynmono,pyrs));
            pvripplepeth(npv,:)=rippeth.ripplepeth.rate(aacs(npv),:);
        end
        ripmodpremono(isnan(ripmodpremono))=0;
        ripmodpremonoCE(isnan(ripmodpremonoCE))=0;
        ripmodDpremonoCE(isnan(ripmodDpremonoCE))=0;
        ripmodSpremonoCE(isnan(ripmodSpremonoCE))=0;
        ripmodpremonoMAX(isnan(ripmodpremonoMAX))=0;
        ripmodpremonoMIN(isnan(ripmodpremonoMIN))=0;
        ripmodpremonoBL(isnan(ripmodpremonoBL))=0;
        ripmodpremonoPETH(isnan(ripmodpremonoPETH))=0;
        ripratepowcorrpremono(isnan(ripratepowcorrpremono))=0;
        ripspikephasemagpremono(isnan(ripspikephasemagpremono))=0;
        ripspikephaseanglepremono(isnan(ripspikephaseanglepremono))=0;
        cumul_numpremono_pv        = [cumul_numpremono_pv, numpremono];
        cumul_numpremono_pvnorm    = [cumul_numpremono_pvnorm, numpremononorm];
        cumul_ripmodpremono_pv     = [cumul_ripmodpremono_pv, ripmodpremono];
        cumul_ripmodpremonoCE_pv   = [cumul_ripmodpremonoCE_pv, ripmodpremonoCE];
        cumul_ripmodNoPre_pv       = [cumul_ripmodNoPre_pv, ripmodNoPre];
        cumul_ripmodDpremonoCE_pv   = [cumul_ripmodDpremonoCE_pv, ripmodDpremonoCE];
        cumul_ripmodSpremonoCE_pv   = [cumul_ripmodSpremonoCE_pv, ripmodSpremonoCE];
        cumul_ripmodpremonoMAX_pv  = [cumul_ripmodpremonoMAX_pv, ripmodpremonoMAX];
        cumul_ripmodpremonoMIN_pv  = [cumul_ripmodpremonoMIN_pv, ripmodpremonoMIN];
        cumul_ripmodpremonoBL_pv   = [cumul_ripmodpremonoBL_pv, ripmodpremonoBL];
        cumul_ripmodpremonoPETH_pv = [cumul_ripmodpremonoPETH_pv, ripmodpremonoPETH];
        cumul_ripmodpremonoMODIFIED = [cumul_ripmodpremonoMODIFIED, ripmodpremonoMODIFIED];
        cumul_ripratepowcorrpremono = [cumul_ripratepowcorrpremono, ripratepowcorrpremono];
        cumul_ripspikephasemagpremono   = [cumul_ripspikephasemagpremono, ripspikephasemagpremono];
        cumul_ripspikephaseanglepremono = [cumul_ripspikephaseanglepremono, ripspikephaseanglepremono];   
        cumul_FRpremono         = [cumul_FRpremono, FRpremono];
        cumul_burstpremono      = [cumul_burstpremono, burstpremono];
        cumul_thetapremono      = [cumul_thetapremono, thetapremono];
        cumul_thetaanglepremono=[cumul_thetaanglepremono thetaanglepremono];
        cumul_thetaanglepremononorm=[cumul_thetaanglepremononorm thetaanglepremononorm];
        cumul_thetamagpremono=[cumul_thetamagpremono thetamagpremono];
        cumul_pvripplepeth.rate     = [cumul_pvripplepeth.rate; pvripplepeth];
        cumul_pvripplepeth.timeEdges= rippeth.ripplepeth.timeEdges;
        cumul_pvThetaHist      =[cumul_pvThetaHist; pvThetaHist];
        cumul_meanpremonotransprob  =[cumul_meanpremonotransprob, meanpremonotransprob];
        cumul_maxpremonotransprob  =[cumul_maxpremonotransprob, maxpremonotransprob];
        
    end

    end

