
% % % % % % % % % % % % % % % % % % % % % % % %
% % % Get cumulative values AACs for
% % % theta phase and ripple modulation
% % % % % % % % % % % % % % % % % % % % % % % %

plotCount = 0;

cumul_ripphase_aac = [];
cumul_ripmod_aac = [];
cumul_ripmodCE_aac      = [];
cumul_ripmodMAX_aac     = [];
cumul_ripmodMIN_aac     = [];
cumul_ripmodBL_aac      = [];
cumul_ripmodPETH_aac    = [];
cumul_ripmod_pyr = [];
cumul_ripmodCE_pyr      = [];
cumul_ripmodMAX_pyr     = [];
cumul_ripmodMIN_pyr     = [];
cumul_ripmodBL_pyr      = [];
cumul_ripmodPETH_pyr    = [];
cumul_ID=[];
cumul_numpremono_aac=[];
cumul_numpremono_aacnorm=[];
cumul_ripmodpremono_aac=[];
cumul_ripmodpremonoCE_aac   = [];
cumul_ripmodpremonoMAX_aac  = [];
cumul_ripmodpremonoMIN_aac  = [];
cumul_ripmodpremonoBL_aac   = [];
cumul_ripmodpremonoPETH_aac = [];
cumul_FRpremono=[];
cumul_thetapremono=[];
cumul_AACripplepeth.rate=[];
cumul_AACripplepeth1ms.rate=[];
cumul_AACThetaHist=[];
cumul_meanpremonotransprob=[];
cumul_maxpremonotransprob=[];
cumul_ripmodpremonoMODIFIED=[];
cumul_FR_aac=[];
IDCt = 0;
MasAACPairSpkC=[];
cumul_RipCCG_aac=[];
cumul_ripRate_pyr = [];
cumul_ripGain_pyr = [];
cumul_ripRate_aac = [];
cumul_ripGain_aac = [];   
cumul_FR_pyr=[];
cumul_ripspkmean_aac    = [];
cumul_ripspkmax_aac     = [];
        
for iSess = sessions
    clear optmod ph_mod STP cell_metrics
    numpremono=[];
    numpremononorm=[];
    ripmodpremonoCE=[];
    ripmodpremono=[];
    ripmodpremonoMAX=[];
    ripmodpremonoMIN=[];
    ripmodpremonoBL=[];
    ripmodpremonoPETH=[];
    ripmodpremonoMODIFIED=[];
    FRpremono=[];
    AACripplepeth=[];
    AACThetaHist=[];
    meanpremonotransprob=[];
    maxpremonotransprob=[];
    AACPairSpkc=[];
    cd(dirN{iSess})
    basepath    = cd;
    basename    = bz_BasenameFromBasepath(basepath);
    load([basename '.ripmod.mat']);
    load([basename '.cell_metrics.cellinfo.mat']);
    load([basename '_celltypes.mat']);
    load([basename '.ripple_ccg.mat'])
    load([basename '.ripspikes.NoStim.analysis.mat'])
    rippeth=load([basename '.ripplepeth.analysis.mat']);
    rippeth1ms=load([basename '.ripplepeth1ms.analysis.mat']);
    mean_ripRate = (mean(ripspikesNoStim.rateperRip_OFF,2))';
    if ~isempty(aacs)
        cumul_ripmodCE_aac      = [cumul_ripmodCE_aac,  cell_metrics.ripples_modulationIndex(aacs)];
        cumul_ripmod_aac        = [cumul_ripmod_aac,  ripmod.mod(aacs)];
        cumul_ripmodMAX_aac     = [cumul_ripmodMAX_aac, ripmod.max(aacs)];
        cumul_ripmodMIN_aac     = [cumul_ripmodMIN_aac, ripmod.min(aacs)];
        cumul_ripmodBL_aac      = [cumul_ripmodBL_aac, ripmod.meanbaseline(aacs)];
        cumul_ripmodPETH_aac    = [cumul_ripmodPETH_aac, ripmod.modPETH(aacs)];
        cumul_FR_aac            = [cumul_FR_aac, cell_metrics.firingRate(aacs)];
        cumul_RipCCG_aac        = [cumul_RipCCG_aac, squeeze(ripple_ccg.ccg(:,end,aacs))];
        cumul_ripRate_aac       = [cumul_ripRate_aac, mean_ripRate(aacs)];
        cumul_ripGain_aac       = [cumul_ripGain_aac, ripspikesNoStim.gainRip_OFF(aacs)];
        ripspikesNoStim.numSpkperRip_OFF(ripspikesNoStim.numSpkperRip_OFF==0)=nan;
        cumul_ripspkmean_aac    = [cumul_ripspkmean_aac nanmean(ripspikesNoStim.numSpkperRip_OFF(aacs,:),2)'];
        cumul_ripspkmax_aac     = [cumul_ripspkmax_aac max(ripspikesNoStim.numSpkperRip_OFF(aacs,:),[],2)']
        for iAAC = aacs;
            IDCt = IDCt +1;
            cumul_ID            = [cumul_ID {[num2str(iSess) '_' num2str(iAAC)]}];
        end
        cumul_ripmodCE_pyr      = [cumul_ripmodCE_pyr,  cell_metrics.ripples_modulationIndex(pyrs)];
        cumul_ripmod_pyr        = [cumul_ripmod_pyr,  ripmod.mod(pyrs)];
        cumul_ripmodMAX_pyr     = [cumul_ripmodMAX_pyr, ripmod.max(pyrs)];
        cumul_ripmodMIN_pyr     = [cumul_ripmodMIN_pyr, ripmod.min(pyrs)];
        cumul_ripmodBL_pyr      = [cumul_ripmodBL_pyr, ripmod.meanbaseline(pyrs)];
        cumul_ripmodPETH_pyr    = [cumul_ripmodPETH_pyr, ripmod.modPETH(pyrs)];
        cumul_ripRate_pyr       = [cumul_ripRate_pyr, mean_ripRate(pyrs)];
        cumul_ripGain_pyr       = [cumul_ripGain_pyr, ripspikesNoStim.gainRip_OFF(pyrs)];
        cumul_FR_pyr            = [cumul_FR_pyr, cell_metrics.firingRate(pyrs)];

        for nAAC = 1:length(aacs);
            transprobs=cell_metrics.putativeConnections.excitatoryTransProb(find(cell_metrics.putativeConnections.excitatory(:,2)==aacs(nAAC)));
            meanpremonotransprob(nAAC)=nanmean(transprobs);
            if isempty(transprobs);
                maxpremonotransprob(nAAC)=NaN;
            else
            maxpremonotransprob(nAAC)=max(transprobs);
            end
            numpremono(nAAC)=length(find(cell_metrics.putativeConnections.excitatory(:,2)==aacs(nAAC)));
            numpremononorm(nAAC)=(length(find(cell_metrics.putativeConnections.excitatory(:,2)==aacs(nAAC)))/length(pyrs));
            presynmono=cell_metrics.putativeConnections.excitatory((find(cell_metrics.putativeConnections.excitatory(:,2) == aacs(nAAC))));
            presynmono=presynmono(ismember(presynmono,pyrs));
            ripmodpremono(nAAC)=nanmean(ripmod.mod(presynmono));
            ripmodpremonoMAX(nAAC)=nanmean(ripmod.max(presynmono));
            ripmodpremonoMIN(nAAC)=nanmean(ripmod.min(presynmono));
            ripmodpremonoBL(nAAC)=nanmean(ripmod.meanbaseline(presynmono));
            ripmodpremonoPETH(nAAC)=nanmean(ripmod.modPETH(presynmono));
            ripmodpremonoMODIFIED(nAAC)=nanmean((ripmod.max(presynmono)-ripmod.min(presynmono))/ripmod.meanbaseline(presynmono));
            FRpremono(nAAC)=mean(cell_metrics.firingRate(presynmono(ismember(presynmono,pyrs))));
            AACripplepeth(nAAC,:)=rippeth.ripplepeth.rate(aacs(nAAC),:);
            AACripplepeth1ms(nAAC,:)=rippeth1ms.ripplepeth.rate(aacs(nAAC),:);
        end
        ripmodpremono(isnan(ripmodpremono))=0;
        ripmodpremonoCE(isnan(ripmodpremonoCE))=0;
        ripmodpremonoMAX(isnan(ripmodpremonoMAX))=0;
        ripmodpremonoMIN(isnan(ripmodpremonoMIN))=0;
        ripmodpremonoBL(isnan(ripmodpremonoBL))=0;
        ripmodpremonoPETH(isnan(ripmodpremonoPETH))=0;
        cumul_numpremono_aac        = [cumul_numpremono_aac, numpremono];
        cumul_numpremono_aacnorm    = [cumul_numpremono_aacnorm, numpremononorm];
        cumul_ripmodpremono_aac     = [cumul_ripmodpremono_aac, ripmodpremono];
        cumul_ripmodpremonoCE_aac   = [cumul_ripmodpremonoCE_aac, ripmodpremonoCE];
        cumul_ripmodpremonoMAX_aac  = [cumul_ripmodpremonoMAX_aac, ripmodpremonoMAX];
        cumul_ripmodpremonoMIN_aac  = [cumul_ripmodpremonoMIN_aac, ripmodpremonoMIN];
        cumul_ripmodpremonoBL_aac   = [cumul_ripmodpremonoBL_aac, ripmodpremonoBL];
        cumul_ripmodpremonoPETH_aac = [cumul_ripmodpremonoPETH_aac, ripmodpremonoPETH];
        cumul_ripmodpremonoMODIFIED = [cumul_ripmodpremonoMODIFIED, ripmodpremonoMODIFIED];
        cumul_FRpremono         = [cumul_FRpremono, FRpremono];
        cumul_AACripplepeth.rate     = [cumul_AACripplepeth.rate; AACripplepeth];
        cumul_AACripplepeth.timeEdges= rippeth.ripplepeth.timeEdges;
        cumul_AACripplepeth1ms.rate     = [cumul_AACripplepeth1ms.rate; AACripplepeth1ms];
        cumul_AACripplepeth1ms.timeEdges= rippeth1ms.ripplepeth.timeEdges;
        cumul_meanpremonotransprob  =[cumul_meanpremonotransprob, meanpremonotransprob];
        cumul_maxpremonotransprob  =[cumul_maxpremonotransprob, maxpremonotransprob];
        
    end

    end

