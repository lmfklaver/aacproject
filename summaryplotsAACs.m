%%<<<<<<< Updated upstream
%% This script makes the summary plots for each AAC

% Create a variable called dirN to loop over to make summary plots:
    launchDirNforAACSessions

    sessions = [1:19];%26:34];%13,14,15,16,17];%1,2,3,4,5];%,8,9,16,17];% 17;%[8,9,16,17]; %1,2,3,4,5,[1,2,3,4,5,6,7,8,9],
%Arch sessions [13,14,15,16,17];%
    %% Get all the gray dots for the theta phase x ripple mod plot via:
    for iSess = sessions
        cd(dirN{iSess})
        basepath = cd;
        basename = bz_BasenameFromBasepath(cd);
        load([basename '.spikes.cellinfo.mat'])
        lfp=bz_GetLFP('all')
        load([basename '_celltypes.mat'])
        load([basename '.gd_eps.mat'])
        load([basename '.ripples.events.mat'])
        [status,interval]=InIntervals(ripples.peaks(:,1),gd_eps); %Detect ripples outside of stim
        ripstart=ripples.timestamps(:,1);
        ripend=ripples.timestamps(:,2);
        gdrips=[];
        gdrips(:,1) = ripstart(status)-.05;
        gdrips(:,2) = ripend(status)+.05;
        [Congdrips] = ConsolidateIntervals(gdrips)
        [SpikeLFPCouplingGdEps]=bz_GenSpikeLFPCoupling(spikes,lfp,'frange',[120 250],'nfreqs',1,'spikeLim',1000000,'cellclass',allcelltypes,'int',Congdrips,'channel',unique([ripples.detectorinfo.detectionparms.channel spikes.maxWaveformCh(aacs)]))
%        [SpikeLFPCouplingStim]=bz_GenSpikeLFPCoupling(spikes,lfp,'frange',[120 250],'spikeLim',8000,'cellclass',allcelltypes,'int',ripspikes.ONrips.timestamps,'channel',unique([ripples.detectorinfo.detectionparms.channel spikes.maxWaveformCh(aacs)]))
        save([basename '.SpikeLFPCouplingGdEps1F.mat'],'SpikeLFPCouplingGdEps')
%        save([basename '.SpikeLFPCouplingStim.mat'],'SpikeLFPCouplingStim')
     end
   for iSess = sessions
        cd(dirN{iSess})
        basepath = cd;
        basename = bz_BasenameFromBasepath(cd);
        load([basename '.gd_eps.mat'])
        load([basename '.spikes.cellinfo.mat'])
        [ripple_ccg] = getRipCCGFixed(basepath,spikes,'epochs',gd_eps,'ccgbin', 0.001,'ccgdur', 1,'saveMat',false);
        save([basename '.ripple_ccg1ms.analysis.mat'], 'ripple_ccg')
   end
  for iSess = sessions
        cd(dirN{iSess})
        basepath = cd;
        basename = bz_BasenameFromBasepath(cd);
        [ripspikesNoStim] = getNumSpkRipNoStim(basepath,'units','all','saveMat',true);
   end
   for iSess = sessions
        cd(dirN{iSess})
        basepath = cd;
        basename = bz_BasenameFromBasepath(cd);
        if isfile([basename '.ripspikes.allripinstim.analysis.mat'])
        load([basename '.ripspikes.allripinstim.analysis.mat'])
        load([basename '.spikes.cellinfo.mat'])
        [zeta] = runZeta(basepath,ripspikes.OFFrips.timestamps(:,1),spikes,'saveMat',true,'saveAs','.ripzeta.stats.mat');
        else
        end
   end
for iSess = sessions
        cd(dirN{iSess})
        basepath = cd;
        basename = bz_BasenameFromBasepath(cd);
        load([basename '.gd_eps.mat'])
        [ripspiketime] = getRipSpkTime(basepath,gd_eps,'units','all','saveMat',true)
        %[ripSTA] = getRipSTA(basepath,gd_eps,'units','all','saveMat',true)
         %if isfile([basename '.ripspikes.allripinstim.analysis.mat'])
        %[ripspiketimeSTIM] = getRipSpkTimeSTIM(basepath,'units','all','saveMat',true)
        %[ripSTASTIM] = getRipSTASTIM(basepath,'units','all','saveMat',true)
         %else
         %end         
end
tic
for iSess = sessions
        cd(dirN{iSess})
        basepath = cd;
        basename = bz_BasenameFromBasepath(cd);
        load([basename '.gd_eps.mat'])
        load([basename '.ripples.events.mat'])
        [status]=InIntervals(ripples.peaks,gd_eps);
        gd_ripplepeaks=ripples.peaks(status);
        [ripplepeth5s] = getPETH_epochs(basepath,'epochs',gd_ripplepeaks,'timwin',[-5 5], ...
                             'binSize', 0.1,'long',true);
        save([basename '.ripplepeth5s.analysis.mat'], 'ripplepeth5s') ;
        toc
end
toc
tic
    for iSess = sessions
        cd(dirN{iSess})
        basepath = cd;
        basename = bz_BasenameFromBasepath(cd);
        load([basename '.ripples.events.mat'])
        load([basename '.optoStim.manipulation.mat'])
        load([basename '.gd_eps.mat'])
        ONrips.detectorinfo=ripples.detectorinfo;
        OFFrips.detectorinfo=ripples.detectorinfo;
        optoStim.detectorinfo=ripples.detectorinfo;
        lfp = bz_GetLFP(ripples.detectorinfo.detectionparms.channel);
        optoStim.peaks = interp1(lfp.timestamps,lfp.timestamps,optoStim.timestamps(:,1),'nearest');
        lfp=[];
        [status]=InIntervals(ripples.peaks,optoStim.timestamps);
        [status1]=InIntervals(ripples.peaks,gd_eps);
        ONrips.peaks=ripples.peaks(status);
        OFFrips.peaks=ripples.peaks(status1);
        if ~exist([basename 'optoStim.lin.evtPSD.mat'])
        getEventPSD(cd,optoStim,ripples.detectorinfo.detectionparms.channel,'space','lin');
        else
        end
%         if sum(status)>0
%         getEventPSD(cd,ONrips,ripples.detectorinfo.detectionparms.channel);
%         getEventPSD(cd,OFFrips,ripples.detectorinfo.detectionparms.channel);
%         else
        %end
        toc
    end
toc

tic
    for iSess = sessions
        cd(dirN{iSess})
        basepath = cd;
        basename = bz_BasenameFromBasepath(cd);
        [thetaEpochs] = detectThetaEpochs('bandpass',[4 10],'powerThreshold',1.5,'force',true);
        load([basename '.ripples.events.mat'])
        load([basename '.optoStim.manipulation.mat'])
        load([basename '.thetaEpochs.states.mat'])
        lfp = bz_GetLFP(thetaEpochs.channel);
        optoStim.peaks = interp1(lfp.timestamps,lfp.timestamps,optoStim.timestamps(:,1),'nearest');
        lfp=[];
        [status]=InIntervals(optoStim.peaks,thetaEpochs.intervals);
        ThetaStim.peaks=optoStim.peaks(status);
        ContStim.peaks=optoStim.peaks(~status);
        getEventPSD(cd,ThetaStim,thetaEpochs.channel);
        getEventPSD(cd,ContStim,thetaEpochs.channel);
        toc
    end
toc

sessions = [24:34];
tic
    for iSess = sessions
        cd(dirN{iSess})
        basepath = cd;
        basename = bz_BasenameFromBasepath(cd);
        load([basename '.thetaEpochs.states.mat'])
        [thetapeth] = getPETH_epochs(basepath,'epochs',thetaEpochs.intervals(:,1),'timwin',[-5 5], ...
                       'binSize', 0.01,'long',true);
        save([basename '.thetapeth.analysis.mat'], 'thetapeth') ;
        toc
    end
toc
tic
for iSess=[13:17]
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(cd);
    load([basename '.ripples.events.mat'])
    load([basename '.optoStim.manipulation.mat'])
    load([basename '_celltypes.mat'])
    load([basename '.rankStats.mat'])
    [eventIDs]=InIntervals(ripples.peaks,optoStim.timestamps);
    figure,subplot(1,3,1),histogram(rankStats.rankClusters(~logical(eventIDs)),'Normalization','Probability')
    title('ID Cluster Stim Ripples')
    subplot(1,3,2)
    histogram(rankStats.rankClusters(logical(eventIDs)),'Normalization','Probability')
    title('ID Cluster Control Ripples')
    subplot(1,3,3)
    plot(rankStats.rankClusters)
    title('SPW-R Cluster ID')
    toc
end
toc
tic
for iSess=sessions
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(cd);
    load([basename '.optoStim.manipulation.mat']);
    computePhaseModulation('excludeIntervals',optoStim.timestamps)
    toc
    close all
end
toc
tic
fracthetarips=[];
for iSess=sessions
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(cd);
    load([basename '.ripples.events.mat'])
    load([basename '.thetaEpochs.states.mat'])
    load([basename '.gd_eps.mat'])
    [status,interval]=InIntervals(ripples.peaks(:,1),gd_eps); %Detect ripples outside of stim
    gd_rips=ripples.peaks(status)
    [events]=InIntervals(gd_rips,thetaEpochs.intervals);
    fracthetarips=[fracthetarips (sum(events)/length(ripples.peaks))]
    toc
    close all
end
toc
tic
numstimrips=[];
numcontrips=[];
for iSess=sessions
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(cd);
    load([basename '.ripples.events.mat'])
    load([basename '.optoStim.manipulation.mat'])
    load([basename '.gd_eps.mat'])
    [status]=InIntervals(ripples.peaks,optoStim.timestamps);
    [status1]=InIntervals(ripples.peaks,gd_eps);
    numstimrips=[numstimrips sum(status)]
    numcontrips=[numcontrips sum(status1)]
end

tic
sessions=[1:19]
distances=[];
aaccount=[];
for iSess=sessions
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(cd);
    load([basename '.cell_metrics.cellinfo.mat'])
    load([basename '_celltypes.mat'])
    [dist,idx]=min(abs(cell_metrics.general.SWR.channelDistance))
    distances=[distances dist]
    aaccount=[aaccount length(aacs)]
    %computePhaseModulation('rippleChannel',ripchan,'plotting',false)
    toc
end
toc


sessions=[1:19]
for iSess=sessions
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(cd);
    load([basename '.spikes.cellinfo.mat'])
    load([basename '_celltypes.mat'])
    load([basename 'ripples.events.mat'])
    for iCell=1:length(spikes.times)
    end
end
toc

% for iSess = sessions
%     cd(dirN{iSess})
%     basepath = cd;
%     basename = bz_BasenameFromBasepath(cd);
%     load([basename '.cell_metrics.cellinfo.mat'])
%     load([basename '.ripple_ccgOFF.mat'])
%     load([basename '_celltypes.mat'])
%     load([basename '.ripmod.mat'])
%     load([basename '.ripple_ccg.mat'])
%     load([basename '.ripple_ccgCont.mat'])
% for iAAC=1:length(aacs)
% presynmono=cell_metrics.putativeConnections.excitatory((find(cell_metrics.putativeConnections.excitatory(:,2) == aacs(iAAC))));
% presynmono=presynmono(ismember(presynmono,pyrs));
% if ~isempty(presynmono)
%     figure,
% subplot(2,2,1)
% hold on
% for monosyn=1:length(presynmono)
% h3 = bar(ripple_ccgOFF.t,ripple_ccgOFF.ccg(:,aacs(iAAC),monosyn));
%             h3.EdgeColor = 'none';
% %            h3.FaceColor = 'r';
%             h3.FaceAlpha = 0.4;
%             title(['CCG within Rip for Presynaptic'])
%             xlabel('time lag(s)')
%             ylabel('Correlation Rate (spikes/second)')
%         %         ylabel('rate')
%         box off
%         %         title('Ripple CCG')
%         legend({'ccg to aac spike'})
%         legend('boxoff')
% end
% [~,n]=find(pyrs==presynmono)
% nonmono=pyrs
% nonmono(n)=[]
% subplot(2,2,3)
% hold on
% for inonmono=1:length(nonmono)
% h3 = bar(ripple_ccgOFF.t,ripple_ccgOFF.ccg(:,aacs(iAAC),nonmono(inonmono)));
%             h3.EdgeColor = 'none';
% %            h3.FaceColor = 'r';
%             h3.FaceAlpha = 0.4;
%             title(['CCG within Rip for nonconnected'])
%             xlabel('time lag(s)')
%             ylabel('Correlation Rate (spikes/second)')
%         %         ylabel('rate')
%         box off
%         %         title('Ripple CCG')
%         legend({'ccg to aac spike'})
%         legend('boxoff')
% end
% subplot(2,2,2)
% hold on
% for monosyn=1:length(presynmono)
% h3 = bar(ripple_ccgCont.t,ripple_ccgCont.ccg(:,monosyn,aacs(iAAC)));
%             h3.EdgeColor = 'none';
% %            h3.FaceColor = 'r';
%             h3.FaceAlpha = 0.4;
%             title(['CCG anytime for Presynaptic'])
%             xlabel('time lag(s)')
%             ylabel('Correlation Rate (spikes/second)')
%             xlim([-.05 .05])
%         %         ylabel('rate')
%         box off
%         %         title('Ripple CCG')
%         legend({'ccg to aac spike'})
%         legend('boxoff')
% end
% [~,n]=find(pyrs==presynmono)
% nonmono=pyrs
% nonmono(n)=[]
% subplot(2,2,4)
% hold on
% for inonmono=1:length(nonmono)
% h3 = bar(ripple_ccgCont.t,ripple_ccgCont.ccg(:,nonmono(inonmono),aacs(iAAC)));
%             h3.EdgeColor = 'none';
% %            h3.FaceColor = 'r';
%             h3.FaceAlpha = 0.4;
%             title(['CCG anytime for nonconnected'])
%             xlabel('time lag(s)')
%             ylabel('Correlation Rate (spikes/second)')
%             xlim([-.05 .05])
%         %         ylabel('rate')
%         box off
%         %         title('Ripple CCG')
%         legend({'ccg to aac spike'})
%         legend('boxoff')
% end
% % suptitle(['Session ' num2str(iSess) 'AAC ' num2str(aacs(iAAC)) 'Ripmod = ' num2str(ripmod.mod(aacs(iAAC)))])
% % unitStr = ['D:\Data\Sorting\CCG_AAC_' num2str(iSess) '_' num2str(iAAC)];
% %             savefig(gcf,[unitStr '.fig'])
% %             print(gcf,[unitStr '.pdf'],'-dpdf','-bestfit')
% %             %         append_pdfs(['E:\Dropbox\PD_Hpc\Progress\AAC\AAC_SummaryPlots_andCCG\SummaryPlot_AACs_20201218_speed2cms.pdf'],[unitStr '.pdf'])
% %             append_pdfs(['D:\Data\Sorting\SummaryPlots_AACs_CCGCent100.pdf'],[unitStr '.pdf'])
% %             
% %             delete([unitStr '.pdf'])
% %             close gcf
% else
% end
% end
% end

getCumulRipModThetaPhase % requires a variable session to work

indnonzero=find(cumul_ripmodpremono_aac>0);
figure,
subplot(3,2,1)
scatter(cumul_RipTTFS_aac(indnonzero),cumul_ripmodCE_aac(indnonzero),'k.')
ylabel('AAC Ripple Modulation')
xlabel('Time to First Spike AAC')
xlim([0 .05])
lsline
ylim([0 3.5])
[r,p,rlo,rup] = corrcoef(cumul_ripmodCE_aac(indnonzero),cumul_RipTTFS_aac(indnonzero))
hold on
text(.02, max(ylim)*0.9, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off

subplot(3,2,2)
scatter(cumul_ripmodpremonoCE_aac(indnonzero),cumul_RipTTFS_aac(indnonzero),'k.')
ylabel('Time to first spike AAC')
xlabel('RipMod of PrePYR')
[r,p,rlo,rup] = corrcoef(cumul_ripmodpremonoCE_aac(indnonzero),cumul_RipTTFS_aac(indnonzero))
ylim([0 .05])
xlim([0 3])
lsline
hold on
text(1.5, max(ylim)*0.9, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off

subplot(3,2,3)
scatter(cumul_numpremono_aac(indnonzero),cumul_RipTTFS_aac(indnonzero),'k.')
ylabel('Time to first spike of AAC')
xlabel('Number of PrePYR')
ylim([0 .05])
xlim([0 9])
lsline
[r,p,rlo,rup] = corrcoef(cumul_numpremono_aac(indnonzero),cumul_RipTTFS_aac(indnonzero))
hold on
text(4, max(ylim)*0.9, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off

subplot(3,2,4)
scatter(cumul_meanpremonotransprob(indnonzero),cumul_RipTTFS_aac(indnonzero),'k.')
ylabel('Time to first spike of AAC')
xlabel('Mean Transmission Probability')
ylim([0 .05])
xlim([0 .1])
lsline
[r,p,rlo,rup] = corrcoef(cumul_meanpremonotransprob(indnonzero),cumul_RipTTFS_aac(indnonzero))
hold on
text(.04, max(ylim)*0.9, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off

subplot(3,2,5)
scatter(cumul_ripphase_aac(indnonzero),cumul_ripmodCE_aac(indnonzero),'k.')
ylabel('Ripmod AAC')
xlabel('Rip Phase Pref (0 peak)')
xlim([-3.5 3.5])
ylim([0 3.5])
[rho pval] = circ_corrcl(cumul_ripphase_aac(indnonzero),cumul_ripmodCE_aac(indnonzero))
hold on
text(.04, max(ylim)*0.9, sprintf('P Value %0.4f', pval), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off

subplot(3,2,6)
scatter(cumul_ripmodCE_aac(indnonzero),cumul_RipSAP_aac(indnonzero),'k.')
ylabel('Average ST around RipPeak')
xlabel('Ripmod AAC')
ylim([0 0.015])
xlim([0 3.5])
lsline
[r,p,rlo,rup] = corrcoef(cumul_RipSAP_aac(indnonzero),cumul_ripmodCE_aac(indnonzero))
hold on
text(2, max(ylim)*0.9, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off

figure,
subplot(3,2,1)
scatter(cumul_ripmodCE_aac(indnonzero),cumul_numpremono_aac(indnonzero),'k.')
ylabel('Number of presynaptic pyrs')
xlabel('AAC Ripple Modulation')
xlim([0 3])
lsline
subplot(3,2,2)
scatter(cumul_ripmodCE_aac(indnonzero),cumul_numpremono_aacnorm(indnonzero),'k.')
ylabel('Percentage presynaptic pyrs')
xlabel('AAC Ripple Modulation')
xlim([0 3])
lsline
subplot(3,2,3)
scatter(cumul_ripmodCE_aac(indnonzero),cumul_FRpremono(indnonzero),'k.')
ylabel('FR of Presynaptic PYR')
xlabel('AAC Ripple Modulation')
xlim([0 3])
lsline
subplot(3,2,4)
scatter(cumul_ripmodCE_aac(indnonzero),cumul_ripmodpremonoCE_aac(indnonzero),'k.')
ylabel('Presynaptic Ripple Modulation')
xlabel('AAC Ripple Modulation')
xlim([0 3])
lsline
subplot(3,2,5)
scatter(cumul_ripmodCE_aac(indnonzero),cumul_burstpremono(indnonzero),'k.')
ylabel('Presynaptic Burstiness')
xlabel('AAC Ripple Modulation')
xlim([0 3])
lsline
subplot(3,2,6)
scatter(cumul_ripmodCE_aac(indnonzero),cumul_thetapremono(indnonzero),'k.')
ylabel('Presynaptic ThetaModIdx')
xlabel('AAC Ripple Modulation')
xlim([0 3])
lsline
scatter(cumul_thetamod_aac(indnonzero),cumul_thetapremono(indnonzero),'k.')
ylabel('Presynaptic ThetaModIdx')
xlabel('AAC Ripple Modulation')
xlim([0 3])
lsline

figure,
scatter(cumul_ripmodpremonoCE_aac(indnonzero),cumul_ripmodCE_aac(indnonzero),'k.')
xlabel('Presynaptic Ripple Modulation')
ylabel('AAC Ripple Modulation')
ylim([0 3.5])
title('AAC Rip Mod x Presynaptic Rip Mod')
xlim([0 3.5])
lsline
[r,p]=corr(cumul_ripmodpremonoCE_aac(indnonzero)',cumul_ripmodCE_aac(indnonzero)','type','Spearman')
hold on
scatter(cumul_ripmodNoPre_aac(indnonzero),cumul_ripmodCE_aac(indnonzero),'r.')
xlim([0 3.5])
lsline
[r,p]=corr(cumul_ripmodNoPre_aac(indnonzero)',cumul_ripmodCE_aac(indnonzero)','type','Spearman')

figure,
scatter(cumul_numpremono_aac(indnonzero),cumul_ripmodCE_aac(indnonzero),'k.')
xlabel('Number of Presynaptic Partners')
ylabel('AAC Ripple Modulation')
ylim([0 3.5])
xlim([0 8.5])
title('AAC Rip Mod x Num Premono')
lsline
[r,p]=corr(cumul_numpremono_aac(indnonzero)',cumul_ripmodCE_aac(indnonzero)','type','Spearman')

figure,
scatter(cumul_FRpremono(indnonzero),cumul_ripmodCE_aac(indnonzero),'k.')
xlabel('Presynaptic Firing Rate (hz)')
ylabel('AAC Ripple Modulation')
ylim([0 3.5])
title('AAC Rip Mod x Presynaptic FR')
xlim([0 3.5])
lsline
[r,p]=corr(cumul_FRpremono(indnonzero)',cumul_ripmodCE_aac(indnonzero)','type','Spearman')

Sindnonzero=find(cumul_ripmodSpremonoCE_aac>0);
figure,
scatter(cumul_ripmodSpremonoCE_aac(Sindnonzero),cumul_ripmodCE_aac(Sindnonzero),'k.')
xlabel('Presynaptic Ripple Modulation (CE)')
ylabel('AAC Ripple Modulation (CE)')
ylim([0 3.5])
title('AAC Rip Mod x Presynaptic Rip Mod')
xlim([0 3.5])
lsline
[r,p]=corrcoef(cumul_ripmodSpremonoCE_aac(Sindnonzero),cumul_ripmodCE_aac(Sindnonzero))

Dindnonzero=find(cumul_ripmodDpremonoCE_aac>0);
figure,
scatter(cumul_ripmodDpremonoCE_aac(Dindnonzero),cumul_ripmodCE_aac(Dindnonzero),'k.')
xlabel('Presynaptic Deep Ripple Modulation (CE)')
ylabel('AAC Ripple Modulation (CE)')
ylim([0 3.5])
title('AAC Rip Mod x Presynaptic Deep Rip Mod')
xlim([0 3.5])
lsline
[r,p]=corrcoef(cumul_ripmodDpremonoCE_aac(Dindnonzero),cumul_ripmodCE_aac(Dindnonzero))

figure,
scatter(cumul_ripratepowcorrpremono(indnonzero),cumul_ratepowercorr_aac(indnonzero),'k.')
xlabel('Presynaptic Ripple Rate Power Correlation')
ylabel('AAC Ripple Rate Power Correlation')
xlim([-.1 .5])
lsline
ylim([-.1 .5])
title('AAC Rip Rate Power Corr x Presynaptic Rip Rate Power Corr')
[r,p]=corrcoef(cumul_ripratepowcorrpremono(indnonzero),cumul_ratepowercorr_aac(indnonzero))

figure,
scatter(cumul_ripspikephasemagpremono(indnonzero),cumul_spikephasemag_aac(indnonzero),'k.')
xlabel('Presynaptic Ripple Spike Phase Magnitude')
ylabel('AAC Ripple Spike Phase Magnitude')
xlim([0 .9])
lsline
ylim([0 .9])
title('AAC Rip SPM x Presynaptic Rip SPM')
[r,p]=corrcoef(cumul_ripspikephasemagpremono(indnonzero),cumul_spikephasemag_aac(indnonzero))

figure,
scatter(cumul_ripspikephaseanglepremono(indnonzero),cumul_spikephaseangle_aac(indnonzero),'k.')
xlabel('Presynaptic Ripple Spike Phase Angle')
ylabel('AAC Ripple Spike Phase Magnitude')
xlim([0 .9])
lsline
ylim([0 .9])
title('AAC Rip SPM x Presynaptic Rip SPM')
[r,p]=corrcoef(cumul_ripspikephaseanglepremono(indnonzero),cumul_spikephaseangle_aac(indnonzero))


figure,
scatter(cumul_maxpremonotransprob(indnonzero),cumul_ripmodCE_aac(indnonzero),'k.')
xlabel('Presynaptic max trans probability')
ylabel('AAC Ripple Modulation')
title('AAC Rip Mod x Presynaptic Max Trans Prob')
lsline
[r,p]=corrcoef(cumul_maxpremonotransprob(indnonzero),cumul_ripmodCE_aac(indnonzero))

figure,
scatter(cumul_meanpremonotransprob(indnonzero),cumul_ripmodCE_aac(indnonzero),'k.')
xlabel('Presynaptic mean trans probability')
ylabel('AAC Ripple Modulation')
title('AAC Rip Mod x Presynaptic mean trans probability')
lsline
[r,p]=corrcoef(cumul_meanpremonotransprob(indnonzero),cumul_ripmodCE_aac(indnonzero))

normripmodpremono=cumul_ripmodpremonoCE_aac(indnonzero).*cumul_meanpremonotransprob(indnonzero)
figure,
scatter(normripmodpremono,cumul_ripmodCE_aac(indnonzero),'k.')
xlabel('Presynaptic Normalized Ripple Modulation (MeanRipmod*MeanTransProb)')
ylabel('AAC Ripple Modulation')
title('AAC Rip Mod x Presynaptic Rip Mod')
lsline
[r,p]=corrcoef(normripmodpremono,cumul_ripmodCE_aac(indnonzero))

figure,
scatter(cumul_FR_aac,cumul_ripmodCE_aac,'k.')
xlabel('AAC Firing Rate (hz)')
ylabel('AAC Ripple Modulation')
title('AAC Rip Mod x AAC Firing Rate')
lsline
[r,p]=corrcoef(cumul_FR_aac,cumul_ripmodCE_aac)

figure,
scatter(cumul_thetaphase_aac,cumul_thetamod_aac,'k.')
xlabel('AAC Pref Theta Phase')
ylabel('AAC Theta Modulation')
title('AAC Pref Theta Phase x AAC Theta Mod')
lsline
[rho pval] = circ_corrcl(cumul_thetaphase_aac, cumul_thetamod_aac)

figure,
histogram(cumul_FR_aac,'BinWidth',.5)
zPeth=zscore(cumul_AACripplepeth.rate,[],2);
[sortzPeth,SI]=sortrows(cumul_ripmodCE_aac','descend');
figure,
imagesc(zPeth(SI,:))
xt = get(gca, 'XTick');
xtnew = linspace(0, max(xt),5);                            
xtlbl = linspace(min(cumul_AACripplepeth.timeEdges), max(cumul_AACripplepeth.timeEdges), numel(xtnew));                  
set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
ylabel('Cell Number')
xlabel('Time (s)')
title('AAC Ripple PETH');
colorbar;
colormap(parula);
clim([-3 6])
hold on
linbot=[50 50]';
lintop=[0 18];
line(linbot,lintop,'Color','white')

thetaout=isoutlier(cumul_thetaphase_aac)
ripmodout=isoutlier(cumul_ripmodCE_aac)
outliers=logical(thetaout+ripmodout)
x(:,1)=(cumul_thetaphase_aac'),x(:,2)=(cumul_ripmodCE_aac')
figure,
scatter(cumul_thetaphase_aac,cumul_ripmodCE_aac,500,'k.')
xlim([-pi pi])
lsline
hold on
ft1 = fittype('a*x^2+b*x+c');
fit1 = fit(x(:,1),x(:,2),ft1);
plot(fit1)
methods(fit1)
xint = linspace(min(-pi),max(pi),100);
CIF = predint(fit1,xint,0.95,'Functional');
CIO = predint(fit1,xint,0.95,'obs');
plot(fit1)
hold on
plot(xint,CIF,':b')
plot(xint,CIO,':g')
xlabel('Preferred Theta Phase')
ylabel('SPW-R Modulation Index')
[rho pval] = circ_corrcl(cumul_thetaphase_aac, cumul_ripmodCE_aac)


[rho pval] = circ_corrcl(cumul_ripphase_aac, cumul_ripmodCE_aac)
x=[]
ripout=isoutlier(cumul_ripphase_aac)
thetamodout=isoutlier(cumul_thetamod_aac)
outliers=logical(ripout+thetamodout)
x(:,1)=(cumul_ripphase_aac(~outliers)'),x(:,2)=(cumul_thetamod_aac(~outliers)')
figure,
scatter(cumul_ripphase_aac(~outliers),cumul_thetamod_aac(~outliers),500,'k.')
xlim([-pi pi])
lsline
hold on
ft1 = fittype('a*x^2+b*x+c');
fit1 = fit(x(:,1),x(:,2),ft1);
plot(fit1)
methods(fit1)
xint = linspace(min(-pi),max(pi),100);
CIF = predint(fit1,xint,0.95,'Functional');
CIO = predint(fit1,xint,0.95,'obs');
plot(fit1)
hold on
plot(xint,CIF,':b')
plot(xint,CIO,':g')
xlabel('Preferred Ripple Phase')
ylabel('Theta Modulation Index')

figure,
ph_bin  = linspace(-pi,pi,32)
histcounts(cumul_thetaphase_pyr,ph_bin)
histogram('BinEdges',ph_bin,'BinCounts',ans)

figure,
subplot(3,2,1)
scatter(cumul_ripmodCE_pyr,cumul_ripmod_pyr,'r.')
[r,p,rlo,rup] = corrcoef(cumul_ripmodCE_pyr,cumul_ripmod_pyr)
hold on
text(.02, max(ylim)*.95, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top','Color','r');
hold on
scatter(cumul_ripmodCE_aac,cumul_ripmod_aac,'k.')
ylabel('cumul ripmod')
xlabel('cumul ripmodCE')
lsline
[r,p,rlo,rup] = corrcoef(cumul_ripmod_aac,cumul_ripmodCE_aac)
hold on
text(.02, max(ylim)*0.9, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off

subplot(3,2,3)
scatter((cumul_ripmodMIN_pyr+cumul_ripmodMAX_pyr)./cumul_ripmodBL_pyr,cumul_ripmod_pyr,'r.')
[r,p,rlo,rup] = corrcoef((cumul_ripmodMIN_pyr+cumul_ripmodMAX_pyr)./cumul_ripmodBL_pyr,cumul_ripmod_pyr)
hold on
text(.02, max(ylim)*.95, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top','Color','r');
scatter((cumul_ripmodMIN_aac+cumul_ripmodMAX_aac)./cumul_ripmodBL_aac,cumul_ripmod_aac,'k.')
ylabel('cumul_ripmod')
xlabel('(cumul_ripmodMIN+cumul_ripmodMAX)/cumul_ripmodBL')
lsline
[r,p,rlo,rup] = corrcoef(cumul_ripmod_aac,(cumul_ripmodMIN_aac+cumul_ripmodMAX_aac)./cumul_ripmodBL_aac)
hold on
text(.02, max(ylim)*0.9, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off

subplot(3,2,5)
scatter((cumul_ripmodMIN_pyr+cumul_ripmodMAX_pyr)./cumul_ripmodBL_pyr,cumul_ripmodCE_pyr,'r.')
[r,p,rlo,rup] = corrcoef((cumul_ripmodMIN_pyr+cumul_ripmodMAX_pyr)./cumul_ripmodBL_pyr,cumul_ripmodCE_pyr)
hold on
text(.02, max(ylim)*.95, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top','Color','r');
scatter((cumul_ripmodMIN_aac+cumul_ripmodMAX_aac)./cumul_ripmodBL_aac,cumul_ripmodCE_aac,'k.')
ylabel('cumul_ripmodCE')
xlabel('(cumul_ripmodMIN+cumul_ripmodMAX)/cumul_ripmodBL')
lsline
[r,p,rlo,rup] = corrcoef(cumul_ripmodCE_aac,(cumul_ripmodMIN_aac+cumul_ripmodMAX_aac)./cumul_ripmodBL_aac)
hold on
text(.02, max(ylim)*0.9, sprintf('P Value %0.4f', p(1,2)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off



%% Select which Subplots you want to plot in the Summary
doWaveform      = true;
doACG           = true;
doPETHPulse     = true;
doPETHRip       = true;

doPhaseMap      = true;
doPhaseMapRip   = false;
doZETA          = true;
doPETHRun       = true;
doFRRun         = true;
doRipTheta      = true;
doLatency       = true;
doRippleLong    = true;
doSave          = true;

%% Start building the Figure
for iSess = sessions  
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(cd);
    % % % % % % % % % % % % %
    % % Get AACs
    % % % % % % % % % % % % %
    load([basename '_celltypes.mat'])
    load([basename '.spikes.cellinfo.mat']) % from buzcode
    load([basename '.mono_res.cellinfo.mat']) % through cellexplorer - unsure if used in current iteration
    
    load([basename '.cell_metrics.cellinfo.mat']) % from cell explorer
    load([basename '.ccginout.mat'])
%    load([basename '.ccginout.analysis.mat']) % through getCCGinout.m
    
    load([basename '.optoStim.manipulation.mat']) % from getPulseEpochs, through a wrapper
    %     load([basename '.dblZeta20_100_100.mat'])
    load([basename '.pethzeta.stats.mat']) %through getZeta.m using ZETA toolbox
    
    ph_mod_5_10Hz = load([basename '.ph_mod_5_10Hz.mat']); %through aacphasemapcode
    ph_mod_39_50Hz = load([basename '.ph_mod_39_50Hz.mat']);
    load([basename '.STP.mat']) % through aacphasemapcode
    load([basename '.ripples.events.mat']) % buzcode
    load([basename '.ripple_ccg.mat'])
    load([basename '.gd_eps.mat'])
    %load([basename '_analogin.mat'])
    
    dotSize = 50;
    % % % % % % % % % % % % %
    % % PulseEpochs
    pulseEpochs = optoStim.timestamps;

    % % % % % % % % % % % % %
    % % runEpochs
    
            if exist([basename '.run.states.mat'],'file')
                load([basename '.run.states.mat'])
                minRunLength = 3;
                minRunSpeed = run.detectorinfo.detectionparms.minRunSpeed;
                selRunEpochs = run.epochs(run.epochs(:,2)-run.epochs(:,1)>=minRunLength,:);
    
            else
            end
basepath = cd;    
    % % % % % % % % % % % % %
    % % Ripple LFP Epochs
    fils  = getAllExtFiles(basepath,'rip',1);
    rip = LoadEvents(fils{1});
    
    % pulls the channel from the ripples and loads the xml file
    rippleChan = str2double(rip.description{1}(regexp(rip.description{1},'[0-9]')));
    %lfp = bz_GetLFP(rippleChan);
    
    
    
    %%
    for iAAC =aacs%1:length(spikes.times)
        
        Fh = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf,'PaperOrientation','Landscape')
        set(gcf,'PaperType','A2')
        
        %  plot panel figure
        %%
        % % % % % % % % % % % % %
        % % Waveform
        % % % % % % % % % % % % %
        if doWaveform
            % CellExplorer does the intan conversion factor already
            subplot(6,4,1)
            xForWave = (1:length(cell_metrics.waveforms.raw{iAAC}))/30000*1000;
            plot(xForWave,cell_metrics.waveforms.raw{iAAC})%
            %         ax1 = gca;
            %         ax1.XAxis.Visible = 'off';
            box off
            ylabel('uV')
            xlabel('time (ms)')
            title(['CE classification ' cell_metrics.putativeCellType(iAAC)])
            
            %             axis off
        end
        %%
        % % % % % % % % % % % % %
        % % Autocorrelation
        % % % % % % % % % % % % %
        
        if doACG
            subplot(6,4,[5 9])
            
            ccg = ccginout.ccgOUT;
            t = ccginout.t;
            
            bar(t,ccg(:,iAAC,iAAC),'k')
            box off
            xlabel('time (s)')
            ylabel('rate')
            set(gca,'TickDir','out')
            
            title(['ACG    CE FR=' num2str(cell_metrics.firingRate((iAAC)))])
        end
        %%
        
        
        %%
        % % % % % % % % % % % % %
        % % PETH Pulse
        % % % % % % % % % % % % %
        if doPETHPulse
            subplot(6,4,2)
            
            if exist([basename '.pulsepeth.analysis.mat'],'file')
                load([basename '.pulsepeth.analysis.mat'])
            else
                
                [peth] = getPETH_epochs(basepath,'epochs',pulseEpochs,...
                    'saveAs', '.pulsepeth.analysis.mat');
            end
            
            ratePulse   = pulsepeth.rate;
            countPulse  = pulsepeth.count;
            timeEdges   = pulsepeth.timeEdges;
            
            % Plot
            h1 = histogram('BinEdges',timeEdges, 'BinCounts',ratePulse(iAAC,:));
            hold on
            box off
            set(gca,'TickDir','out')
            title('PETH to Pulse')
            h1.EdgeColor = 'none';
            h1.FaceColor = 'k';
            xlabel('time(s)')
            ylabel('spikes/s')
            xlim(pulsepeth.timwin)
            
            %%
            % % % % % % % % % % % % %
            % % Raster Pulse
            % % % % % % % % % % % % %
            subplot(6,4,[6 10])
            plotSpkOffset = 0;
            
            selTrialsPulse = pulsepeth.trials{iAAC};
            
            for iPulse = 1:length(selTrialsPulse)
                selPulseTr = selTrialsPulse{iPulse};
                plot(selPulseTr',repmat(plotSpkOffset,1,length(selPulseTr)),'k.');
                hold on
                plotSpkOffset = plotSpkOffset+1;
            end
            
            box off
            set(gca,'ydir','reverse')
            ylimits = get(gca,'YLim');
            xlabel('time (s)')
            ylabel('trials')
            %             set(gca,'TickDir','out')
            ylim([ylimits(1) plotSpkOffset])
            xlim(pulsepeth.timwin)
            
        end
        %%
        % % % % % % % % % % % % %
        % % PETH Ripple
        % % % % % % % % % % % % %
        
        if doPETHRip
            
            subplot(6,4,3)
            
            load([basename '.ripples.events.mat'])
            if exist([basename '.ripplepeth.analysis.mat'],'file')
                load([basename '.ripplepeth.analysis.mat'])
            else
            [status]=InIntervals(ripples.peaks,gd_eps);
            gd_ripplepeaks=ripples.peaks(status);
%                 [peth] = getPETH_epochs(basepath,'epochs',gd_ripplepeaks,'timwin',[-0.4 0.4], ...
%                     'binSize', 0.01, 'saveAs', '.ripplepeth.analysis.mat');
                   [peth] = getPETH_epochs(basepath,'epochs',gd_ripplepeaks,'timwin',[-0.5 .5], ...
                    'binSize', 0.01, 'saveAs', '.ripplepeth.analysis.mat');
            end
 
            rateHistoRip    = ripplepeth.rate;
            countPulse      = ripplepeth.count;
            timeEdges       = ripplepeth.timeEdges;
            
            %
            % Plot
            h1 = histogram('BinEdges',timeEdges, ...
                'BinCounts',rateHistoRip(iAAC,:));
            hold on
            box off
            title('PETH Ripple')
            xlabel('time(s)')
            ylabel('spikes/s')
            h1.EdgeColor = 'none';
            h1.FaceColor = 'k';
            xlim(ripplepeth.timwin)
            
            
            %%
            % % % % % % % % % % % % %
            % % Raster Ripple
            % % % % % % % % % % % % %
            
            subplot(6,4,[7 11])
            plotSpkOffset = 0;
            selRipTr=[];
            selTrialsRip = ripplepeth.trials{iAAC};
            [status]=InIntervals(ripples.peaks,gd_eps);
            gd_ripplepeaks=ripples.peaks(status);
            for iRip = 1:length(selTrialsRip)
                selRipTr = selTrialsRip{iRip};
                plot(selRipTr',repmat(plotSpkOffset,1,length(selRipTr)),'k.');
                hold on
                plotSpkOffset = plotSpkOffset+1;
            end
            
            box off
            set(gca,'ydir','reverse')
            ylimits = get(gca,'YLim');
            xlabel('time (s)')
            ylabel('trials')
            set(gca,'TickDir','out')
            ylim([ylimits(1) plotSpkOffset])
            xlim(ripplepeth.timwin)
        end
        %%
        % % % % % % % % % % % % %
        % % PETH RUN Onset
        % % % % % % % % % % % % %
        if doPETHRun
            if ~isempty(selRunEpochs)
                
               
                    subplot(6,4,4)
                    
                    %                     if exist([basename '.runpeth.analysis.mat'],'file')
                    %                         load([basename '.runpeth.analysis.mat'])
                    %                     else
                    
                    [peth] = getPETH_epochs(basepath,'epochs',selRunEpochs,...
                        'timwin',[-5 5],'binSize',0.1,'saveAs','.runpeth2cm.analysis.mat');
                    %                     end
                    
                    rateRun   = peth.rate;
                    countRun  = peth.count;
                    timeEdges   = peth.timeEdges;
                    
                    % Plot
                    h1 = histogram('BinEdges',timeEdges, 'BinCounts',rateRun(iAAC,:));
                    hold on
                    box off
                    set(gca,'TickDir','out')
                    title('PETH to Run Onset')
                    
                    xlabel('time(s)')
                    ylabel('spikes/s')
                    xlim(peth.timwin)
                    h1.EdgeColor = 'none';
                    h1.FaceColor = 'k';
                
                
                
                
                %%
                % % % % % % % % % % % % %
                % % Raster RUN Onset
                % % % % % % % % % % % % %
                
                subplot(6,4,[8,12])
                
                selTrialsRun = peth.trials{iAAC};
                plotSpkOffset = 0;
                
                for iRun = 1:size(selRunEpochs,1)
                    selRunTr = selTrialsRun{iRun};
                    plot(selRunTr',repmat(plotSpkOffset,1,length(selRunTr)),'k.');
                    hold on
                    plotSpkOffset = plotSpkOffset+1;
                end
                
                
                box off
                set(gca,'ydir','reverse')
                ylimits = get(gca,'YLim');
                xlabel('time (s)')
                ylabel('trials')
                %             set(gca,'TickDir','out')
                ylim([ylimits(1) plotSpkOffset])
                xlim(peth.timwin)
            end
        end
        
        
        %%
        % % % % % % % % % % % % %
        % % Monosynaptic Connections
        % % % % % % % % % % % % %
        subplot(6,4,13)
        load([basename '.STP.mat'])
        aacIDX = find(cell_metrics.putativeConnections.excitatory(:,2) == iAAC);
        aacID=(cell_metrics.putativeConnections.excitatory(aacIDX,1));
        plot(cell_metrics.acg.wide(400:600,aacID));
        box off
        xlim([1 201]);
        xIndVals = 50:100:450;
        xLabelSec = xIndVals-((200)/2);% align around center;
        xLabelSec = xLabelSec*1; % (mono_res.binsize)
        set(gca,'XTick',xIndVals,'XTickLabel', num2cell(xLabelSec))
        
        xlabel('time (s)')
        ylabel('rate')
        title('Presynaptic Partners')
        % % % % % % % % % % % % %
        % % Non-Monosynaptic Connections VS monosynaptic for the session
        % % % % % % % % % % % % %
        subplot(6,4,16)
        allpreIDX=[];
        allaacIDX=[];
        for iallaacs=1:length(aacs)
        allpreIDX = find(cell_metrics.putativeConnections.excitatory(:,2) == aacs(iallaacs));
        allaacIDX=[allaacIDX;cell_metrics.putativeConnections.excitatory(allpreIDX,1)];
        end
        not_mono_con_to_aac=unique([cell_metrics.putativeConnections.excitatory(:,1); pyrs'])
        not_mono_con_to_aac_acg=setdiff(not_mono_con_to_aac,allaacIDX);
        uniqueallaacIDX=unique(allaacIDX);
        if ~isempty(uniqueallaacIDX)
        if size(uniqueallaacIDX,1)>1
           meanpresyn=mean(cell_metrics.acg.wide(:,uniqueallaacIDX),2);
        else
            meanpresyn=cell_metrics.acg.wide(:,uniqueallaacIDX);
        end
        else
            continue
        end
        meannopresyn=mean(cell_metrics.acg.wide(:,not_mono_con_to_aac_acg),2);
        plot(meanpresyn(400:600));
        hold on
        plot(meannopresyn(400:600));
        box off
        xlim([1 201])
        xIndVals = 50:100:450;
        xLabelSec = xIndVals-((200)/2);% align around center;
        xLabelSec = xLabelSec*1; % (mono_res.binsize)
        set(gca,'XTick',xIndVals,'XTickLabel', num2cell(xLabelSec))
        
        xlabel('time (s)')
        ylabel('rate')
        title('Connected v Non-connected partners')
        
        %%
        % % % % % % % % % % % % %
        % % CCG Heatmap
        % % % % % % % % % % % % %
        subplot(6,4,17)
        %         load([basename '.ccg.mat'])
        monoIdx = find(cell_metrics.putativeConnections.excitatory(:,2)==iAAC);
        selCCGin = cell_metrics.putativeConnections.excitatory(monoIdx,1);
        selCCGout = iAAC;
        h1 = imagesc(zscore(ccg(:,selCCGin,selCCGout)',[],2));
        %     h1 = imagesc((ccg(:,selCCGin,selCCGout)'));
        
        
        xIndVals = 1:100:201;% 401
        xlim([xIndVals(1) xIndVals(end)]);
        xLabelSec = xIndVals-101;% align around center;201
        xLabelSec = xLabelSec*0.001; % (ccgbinsize ripple)
        set(gca,'XTick',xIndVals,'XTickLabel', num2cell(xLabelSec))
        
        xlabel('time (s)')
        ylabel('Presynaptic Cell')
        cb1= colorbar;
        ylabel(cb1,'Z-scored Rate')
        box off
        
        title('CCG')
        
%         %
%         % % % % % % % % % % % %
%         % ZETA
%         % % % % % % % % % % % %
%                 subplot(6,4,9)
%                 if doZETA
%                 if ~isempty(regexp(basename,'mouse', 'once'))
%                     ZetaP20 = dblZetaPChR20;
%                     ZetaP100 = dblZetaPChR100;
%         
%                 elseif isempty(regexp(basename,'mouse', 'once'))
%                     ZetaP20 = dblZetaPArch20;
%                     ZetaP100 = dblZetaPArch100;
%                 end
%         
%                 dotSize = 50;
%         
%                 scatter(ZetaP20,ZetaP100,dotSize,[211/255,211/255,211/255],'filled')
%                 hold on
%                 scatter(ZetaP20(iAAC),ZetaP100(iAAC),dotSize,'filled','m')
%         
%         
%                 xlabel('p over 20ms');
%                 ylabel('p over 100ms');
%                 xlim([0 1])
%                 ylim([0 1])
%         
%                 xline(0.05,':')
%                 yline(0.05,':')
%         
%         
%                 legend({'all neurons sess','selected AAC'},'Location','northeast','NumColumns',1);
%                 end
%         
%         
        %%
        % % % % % % % % % % % % %
        % % Gain Ripple
        % % % % % % % % % % % % %
        subplot(6,4,15)
        
        h3 = bar(ripple_ccg.t,ripple_ccg.ccg(:,end,iAAC));
            h3.EdgeColor = 'none';
            h3.FaceColor = 'k';
            xlim(ripplepeth.timwin)
            title(['Ripple CCG' num2str(iAAC)])
            xlabel('time lag(s)')
            ylabel('Correlation Rate (spikes/second')
        %         ylabel('rate')
        box off
        %         title('Ripple CCG')
        legend({'ccg to ripple peak'})
        legend('boxoff')
        
        
        
        %%
        % % % % % % % % % % % %
        % Histogram Pref Theta Phase for each Spike
        % % % % % % % % % % % %
        subplot(6,4,23)
        spkInstPhaseStruct=[];
        [spkInstPhaseStruct]=spkInstPhaseOutPulse(basepath, iAAC);
        for q=1:size(spkInstPhaseStruct.binnedHisto,1)
            hold on
            histogram('BinEdges',spkInstPhaseStruct.ph_bin, 'BinCounts',...
                spkInstPhaseStruct.binnedHisto(q,([round(length(spkInstPhaseStruct.binnedHisto)/2)+1:length(spkInstPhaseStruct.binnedHisto) ...
                1:round(length(spkInstPhaseStruct.binnedHisto)/2)])),'DisplayStyle','stairs')
%         legendselfreq{q}= num2str(spkInstPhaseStruct.freqband(q,:));
%         legend({legendselfreq});
%         legend('boxoff')
        end
        histogram('BinEdges',spkInstPhaseStruct.ph_bin, 'BinCounts',spkInstPhaseStruct.meanbinnedHisto([round(length(spkInstPhaseStruct.meanbinnedHisto)/2)+1:length(spkInstPhaseStruct.meanbinnedHisto) ...
                1:round(length(spkInstPhaseStruct.meanbinnedHisto)/2)]));
        hold off
        %%
        % % % % % % % % % % % % %
        % % Ripple Mod x Theta Phase
        % % % % % % % % % % % % %
        if doRipTheta
            subplot(6,4,19)
            hold off
            scatter(cumul_ripphase_aac,  cumul_ripmod_aac,dotSize,[211/255,211/255,211/255],'filled')
            hold on
            
            for iL = 1:length(cumul_ID)
                if find(strcmpi(cumul_ID{iL},[num2str(iSess) '_' num2str(iAAC)]))
                    scatter(cumul_rip_aac(iL),  cumul_ripmod_aac(iL),dotSize,'filled','m')
                end
            end
            xlabel('ripple phase')
            xlim([-pi pi])
            ylabel('ripple mod')
            ylim([0 5])
%             yline(1,':')
            
            lgd = legend({'all AACs','selected AAC'},'Location','northoutside','NumColumns',2);
            legend('boxoff')
        end
        %%
        %%
        % % % % % % % % % % % % %
        % % Latency
        % % % % % % % % % % % % %
        
        if doLatency
            subplot(6,4,14)
            LatencyFirstSpike = zeros(1,length(selTrialsPulse));
            for iPulse = 1:length(selTrialsPulse)
                positiveSpk = selTrialsPulse{iPulse}(selTrialsPulse{iPulse}>0);
                if ~isempty(positiveSpk)
                    LatencyFirstSpike(iPulse) =positiveSpk(1);
                else
                    LatencyFirstSpike(iPulse) =NaN;
                end
            end
            
            edgesLat = 0:0.01:0.3; % NB Hardcoded
            
            CtsLat = histcounts(LatencyFirstSpike, edgesLat);
            histogram('BinEdges', edgesLat, 'BinCounts', CtsLat)
            set(gca,'YScale','log')
            ylabel('count')
            xlabel('latency first spike')
            box off
            
        end

        %%
        
        % % % % % % % % % % % % %
        % % Phasemap
        % % % % % % % % % % % % %
        
        if doPhaseMap
           
            subplot(6,4,18)
            load([basename '.ph_portrait.analysis.mat'])
            ph_bin = linspace(-pi,pi,16);
            k = gaussian2Dfilter([10 10],[.5 .5]);
            
            
            imagesc(ph_portrait.ph_bin,[],nanconvn((ph_portrait.ph_rate(:,[end/2+1:end-1 1:(end/2)],iAAC)),k),...
                    [min(linearize(ph_portrait.ph_rate(:,1:end-1,iAAC))) max(linearize(ph_portrait.ph_rate(:,1:end-1,iAAC)))])            
            hold on
            colormap('jet')
            phasecurv=10+cos(ph_portrait.ph_bin)*10
            plot(ph_portrait.ph_bin(1:end-1),phasecurv([(end/2)+1:end-1 1:end/2]),'w')
             set(gca,'ytick',1:10:length(ph_portrait.freq),'yticklabel',round(ph_portrait.freq(1:10:end)))
            set(gca,'ydir','normal')
%             ylabel('Frequency (logscale)')
            xlabel('Phase')
            title('Phase Portrait')
        end

        
        % % % % % % % % % % % % %
        % % Phasemap Ripple
        % % % % % % % % % % % % %
        
        if doPhaseMapRip
           
            subplot(6,4,18)
            load([basename '.ph_portrait_rip.analysis.mat'])
            ph_bin = linspace(-pi,pi,16);
            k = gaussian2Dfilter([10 10],[.5 .5]);
            
            
            imagesc(ph_portrait.ph_bin,[],nanconvn((ph_portrait.ph_rate(:,[end/2+1:end-1 1:(end/2)],iAAC)),k),...
                    [min(linearize(ph_portrait.ph_rate(:,1:end-1,iAAC))) max(linearize(ph_portrait.ph_rate(:,1:end-1,iAAC)))])            
            hold on
            colormap('jet')
            phasecurv=10+cos(ph_portrait.ph_bin)*10
            plot(ph_portrait.ph_bin(1:end-1),phasecurv([(end/2)+1:end-1 1:end/2]),'w')
             set(gca,'ytick',1:5:length(ph_portrait.freq),'yticklabel',round(ph_portrait.freq(1:5:end)))
            set(gca,'ydir','normal')
            ylabel('Frequency')
            xlabel('Phase')
        title('Phase Portrait')
        
                %%
        % % % % % % % % % % % %
        % Histogram Pref Ripple Phase for each Spike
        % % % % % % % % % % % %
        subplot(6,4,23)
        spkInstPhaseStruct=[];
        [spkInstPhaseStruct]=spkInstPhaseOutPulseRIP(basepath, iAAC)
        for q=1:size(spkInstPhaseStruct.binnedHisto,1)
            hold on
            histogram('BinEdges',spkInstPhaseStruct.ph_bin, 'BinCounts',...
                spkInstPhaseStruct.binnedHisto(q,([round(size(spkInstPhaseStruct.binnedHisto,2)/2)+1:size(spkInstPhaseStruct.binnedHisto,2) ...
                1:round(size(spkInstPhaseStruct.binnedHisto,2)/2)])),'DisplayStyle','stairs')
%         legendselfreq{q}= num2str(spkInstPhaseStruct.freqband(q,:));
%         legend({legendselfreq});
%         legend('boxoff')
        end
        histogram('BinEdges',spkInstPhaseStruct.ph_bin, 'BinCounts',spkInstPhaseStruct.meanbinnedHisto([round(length(spkInstPhaseStruct.meanbinnedHisto)/2)+1:length(spkInstPhaseStruct.meanbinnedHisto) ...
                1:round(length(spkInstPhaseStruct.meanbinnedHisto)/2)]));
        hold off
        end
        % % % % % % % % % % % %
        % % Ripple mod x thetaphase
        % % % % % % % % % % % % %

        
            if doRipTheta
            subplot(6,4,22)
            hold off
            scatter(cumul_thetamod_aac,  cumul_ripmod_aac,dotSize,[211/255,211/255,211/255],'filled')
            hold on
            
            for iL = 1:length(cumul_ID)
                if find(strcmpi(cumul_ID{iL},[num2str(iSess) '_' num2str(iAAC)]))
                    scatter(cumul_thetamod_aac(iL),  cumul_ripmod_aac(iL),dotSize,'filled','m')
                end
            end
            xlabel('theta mod')
            ylabel('ripple mod')
            ylim([0 5])
%             yline(1,':')
            
            lgd = legend({'all AACs','selected AAC'},'Location','northoutside','NumColumns',2);
            legend('boxoff')
        end
        
        %%
        % % % % % % % % % % % % %
        % %Spikes per Ripple Dist
        % % % % % % % % % % % % %
        
        subplot(6,4,21)
        if exist([basename '.ripspikes.allripinstim.analysis.mat'])
            load([basename '.ripspikes.allripinstim.analysis.mat'])
        maxHistoRip = max(ripspikes.spikesRipNum{iAAC});
        edgesRip = 0:maxHistoRip;
        xtRip = [0.5:5:maxHistoRip+0.5];
        histogram(ripspikes.spikesRipNum{iAAC},edgesRip)
        xlabel('Max Num Spk in Rip')
        ylabel('count')
        set(gca,'YScale','log','XTick',xtRip-0.5,'XTickLabel', num2cell(xtRip-0.5))
        box off
        else

        
%         %%
        % % % % % % % % % % % % %
        % %Ripple cycle spikes
        % % % % % % % % % % % % %
        subplot(6,4,23)
        
        
        edgesRip = 0:1:5;
        histogram(cell2mat(ripspikes.numSpkPerCycPerRip(iAAC)),edgesRip)
        xlabel('Avg Number of Spikes per Ripple Cycle')
        ylabel('count')
        set(gca,'YScale','log','XTick',edgesRip,'XTickLabel', num2cell(edgesRip))
        box off
                end
        %%
        %%
        % % % % % % % % % % % % %
        % % Location
        % % % % % % % % % % % % %
        subplot(6,4,[20 24])
        
        
        load('chanMap.mat')
        aacChanMap = chanMap == spikes.maxWaveformCh(iAAC)+1;
        ripChanMap = chanMap == rippleChan+1;
        
        plot(xcoords(ripChanMap)-1,ycoords(ripChanMap),'ko','MarkerFaceColor','k')
        hold on
        plot(xcoords(aacChanMap)-.5,ycoords(aacChanMap),'mo','MarkerFaceColor','m')
        text(xcoords,ycoords,num2str(chanMap0ind))
        xlim([min(xcoords)-10, max(xcoords)+10])
        ylim([min(ycoords)-10 ,max(ycoords)+10])
        
        box off
        xlabel('distance (um)')
        ylabel('distance (um)')
        lloc = legend({'Ripple','AAC'},'Location','northoutside','Box','off','NumColumns',2);
        
        %%
        
        % % % % % % % % % % % % %
        % % Suptitle + Save
        % % % % % % % % % % % % %
        
        suptitle({['Unit ' num2str(iSess) '_' num2str(spikes.UID(iAAC))]})
        
        if doSave
            unitStr = ['D:\Data\Sorting\SummaryPlot_AAC_rip' num2str(iSess) '_' num2str(iAAC)];
            savefig(gcf,[unitStr '.fig'])
            print(gcf,[unitStr '.pdf'],'-dpdf','-bestfit')
            %         append_pdfs(['E:\Dropbox\PD_Hpc\Progress\AAC\AAC_SummaryPlots_andCCG\SummaryPlot_AACs_20201218_speed2cms.pdf'],[unitStr '.pdf'])
            append_pdfs(['D:\Data\Sorting\SummaryPlots_AACs_Newest.pdf'],[unitStr '.pdf'])
            
            delete([unitStr '.pdf'])
            close gcf
        end
        
    end    
    end