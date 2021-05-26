%% This script makes the summary plots for each AAC

% Create a variable called dirN to loop over to make summary plots:
launchDirNforAACSessions

sessions = [1,2,3,4,5,8,9,16,17];% 17;%[8,9,16,17]; %1,2,3,4,5,

% Sess 6 has pulse artifacts
% Sess 7 needs to be checked, was copied over with an empty .dat initially
% Sess 10, 11 and 12 have no good epochs %
% Sess 13 has 200ms pulses
% Sess 14 has no AACs

%% Get all the gray dots for the theta phase x ripple mod plot via:
getCumulRipModThetaPhase % requires a variable session to work


%% Select which Subplots you want to plot in the Summary
doWaveform      = true;
doACG           = true;
doPETHPulse     = true;
doPETHRip       = true;
doPhaseMap      = true;
doZETA          = false;
doPETHRun       = true;
doFRRun         = false;
doRipTheta      = false;
doLatency       = true;
doRippleLong    = false;
doSave          = false;

%% Start building the Figure

for iSess = 1%sessions
    
    cd(dirN{iSess})
    basepath = cd;
    basename = bz_BasenameFromBasepath(cd);
    
    
    % % % % % % % % % % % % %
    % % Get AACs
    % % % % % % % % % % % % %
    
    [~, ~, aacs] = splitCellTypes(basepath);
    
    load([basename '.spikes.cellinfo.mat']) % from buzcode
    load([basename '.mono_res.cellinfo.mat']) % through cellexplorer - unsure if used in current iteration
    
    load([basename '.cell_metrics.cellinfo.mat']) % from cell explorer
    %     load([basename '.ccginout.mat'])
    load([basename '.ccginout.analysis.mat']) % through getCCGinout.m
    
    load([basename '.optoStim.manipulation.mat']) % from getPulseEpochs, through a wrapper
    %     load([basename '.dblZeta20_100_100.mat'])
    load([basename '.pethzeta.stats.mat']) %through getZeta.m using ZETA toolbox
    
    load([basename '.ph_mod.mat']) %through aacphasemapcode
    load([basename '.STP.mat']) % through aacphasemapcode
    load([basename '.ripples.events.mat']) % buzcode
    
    
    
    % % % % % % % % % % % % %
    % % PulseEpochs
    pulseEpochs = optoStim.timestamps;
    
    % % % % % % % % % % % % %
    % % runEpochs
    
    minRunLength = 3;
    minRunSpeed = 2;
    %
    %         if exist([basename '.run.states.mat'],'file')
    %             load([basename '.run.states.mat'])
    %             minRunSpeed = run.detectorinfo.detectionparms.minRunSpeed
    %             selRunEpochs = run.epochs(run.epochs(:,2)-run.epochs(:,1)>=minRunLength,:);
    %
    %         else
    load([basename '_analogin.mat'])
    selRunEpochs = [];
    
    if isfield(analogin,'pos')
        if doPETHRun
            %             if exist([basename '.run2cm.states.mat'], 'file')
            %                 load ([basename '.run2cm.states.mat'])
            %             else
            
            [vel] = getVelocity(analogin,'doFigure',false,'downsampleFactor',3000);
            [run] = getRunEpochs(basepath,vel,'minRunSpeed',minRunSpeed,'saveMat',true,'saveAs','.run2cm.states.mat');
            %             end
            
            selRunEpochs = run.epochs(run.epochs(:,2)-run.epochs(:,1)>=minRunLength,:);
        end
    end
    %     end
    
    % % % % % % % % % % % % %
    % % Ripple LFP Epochs
    fils  = getAllExtFiles(basepath,'rip',1);
    rip = LoadEvents(fils{1});
    
    % pulls the channel from the ripples and loads the xml file
    rippleChan = str2double(rip.description{1}(regexp(rip.description{1},'[0-9]')));
    lfp = bz_GetLFP(rippleChan);
    
    
    
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
            title(['CE classification ' cell_metrics.putativeCellType{(iAAC)}]) %', CE_FR = ' num2str(spikes.CE_FR((iUnit)))])
            
            %             axis off
        end
        %%
        % % % % % % % % % % % % %
        % % Autocorrelation
        % % % % % % % % % % % % %
        
        if doACG
            subplot(6,4,5)
            
            ccg = ccginout.ccgOUT;
            t = ccginout.t;
            
            bar(t,ccg(:,iAAC,iAAC),'k')
            box off
            xlabel('time (s)')
            ylabel('rate')
            set(gca,'TickDir','out')
            
            title('ACG')
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
            
            ratePulse   = peth.rate;
            countPulse  = peth.count;
            timeEdges   = peth.timeEdges;
            
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
            xlim(peth.timwin)
            
            %%
            % % % % % % % % % % % % %
            % % Raster Pulse
            % % % % % % % % % % % % %
            subplot(6,4,[6 10])
            plotSpkOffset = 0;
            
            selTrialsPulse = peth.trials{iAAC};
            
            for iPulse = 1:length(pulseEpochs)
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
            xlim(peth.timwin)
            
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
                [peth] = getPETH_epochs(basepath,'epochs',ripples.peaks,'timwin',[-0.4 0.4], ...
                    'binSize', 0.01, 'saveAs', '.ripplepeth.analysis.mat');
            end
            
            rateHistoRip    = peth.rate;
            countPulse      = peth.count;
            timeEdges       = peth.timeEdges;
            
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
            xlim(peth.timwin)
            
            
            %%
            % % % % % % % % % % % % %
            % % Raster Ripple
            % % % % % % % % % % % % %
            
            subplot(6,4,[7 11])
            plotSpkOffset = 0;
            
            selTrialsRip = peth.trials{iAAC};
            for iRip = 1:length(ripples.timestamps)
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
            xlim(peth.timwin)
        end
        %%
        % % % % % % % % % % % % %
        % % PETH RUN Onset
        % % % % % % % % % % % % %
        if doPETHRun
            if ~isempty(selRunEpochs)
                
                if isfield(analogin,'pos')
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
                end
                
                
                
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
        aacIDX = find(STP.post_idx(:,3) == iAAC);
        plot(STP.acg_pre(aacIDX,:)')
        box off
        xlim([1 size(STP.acg_pre,2)])
        xIndVals = 50:100:450;
        xLabelSec = xIndVals-((size(STP.acg_pre,2)-1)/2);% align around center;
        xLabelSec = xLabelSec*mono_res.binSize; % (mono_res.binsize)
        set(gca,'XTick',xIndVals,'XTickLabel', num2cell(xLabelSec))
        
        xlabel('time (s)')
        ylabel('rate')
        title('Presynaptic Partners')
        
        
        %%
        % % % % % % % % % % % % %
        % % CCG Heatmap
        % % % % % % % % % % % % %
        subplot(6,4,17)
        %         load([basename '.ccg.mat'])
        monoIdx = find(STP.post_idx(:,3)==iAAC);
        selCCGin = STP.pre_idx(monoIdx,3);
        selCCGout = iAAC;
        h1 = imagesc(zscore(ccg(:,selCCGin,selCCGout)',[],2));
        %     h1 = imagesc((ccg(:,selCCGin,selCCGout)'));
        
        
        xIndVals = 1:100:201;% 401
        xlim([xIndVals(1) xIndVals(end)])
        xLabelSec = xIndVals-101;% align around center;201
        xLabelSec = xLabelSec*0.001; % (ccgbinsize ripple)
        set(gca,'XTick',xIndVals,'XTickLabel', num2cell(xLabelSec))
        
        xlabel('time (s)')
        ylabel('Presynaptic Cell')
        cb1= colorbar;
        ylabel(cb1,'Z-scored Rate')
        box off
        
        title('CCG')
        
        %%
        % % % % % % % % % % % % %
        % % ZETA
        % % % % % % % % % % % % %
        %         subplot(6,4,9)
        %         if doZETA
        %         if ~isempty(regexp(basename,'mouse', 'once'))
        %             ZetaP20 = dblZetaPChR20;
        %             ZetaP100 = dblZetaPChR100;
        %
        %         elseif isempty(regexp(basename,'mouse', 'once'))
        %             ZetaP20 = dblZetaPArch20;
        %             ZetaP100 = dblZetaPArch100;
        %         end
        %
        %         dotSize = 50;
        %
        %         scatter(ZetaP20,ZetaP100,dotSize,[211/255,211/255,211/255],'filled')
        %         hold on
        %         scatter(ZetaP20(iAAC),ZetaP100(iAAC),dotSize,'filled','m')
        %
        %
        %         xlabel('p over 20ms');
        %         ylabel('p over 100ms');
        %         xlim([0 1])
        %         ylim([0 1])
        %
        %         xline(0.05,':')
        %         yline(0.05,':')
        %
        %
        %         legend({'all neurons sess','selected AAC'},'Location','northeast','NumColumns',1);
        %         end
        
        
        %%
        % % % % % % % % % % % % %
        % % Gain Ripple
        % % % % % % % % % % % % %
        subplot(6,4,15)
        
        plot(ph_mod.ripple_ccg.ccg(iAAC,:))
        xlim([4900 5100])
        currX = get(gca,'xtick');
        xIndVals = [4900 5000 5100];
        xLabelSec = xIndVals-5000;% align around center;
        xLabelSec = xLabelSec*0.005; % (ccgbinsize ripple)
        set(gca,'XTick',xIndVals,'XTickLabel', num2cell(xLabelSec))
        set(gca,'FontSize',8.5)
        
        xlabel('time (s)')
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
        
        %%
        % % % % % % % % % % % % %
        % % Ripple Mod x Theta Phase
        % % % % % % % % % % % % %
        if doRipTheta
            subplot(6,4,19)
            hold off
            scatter(cumul_thetaphase_aac,  cumul_ripmod_aac,dotSize,[211/255,211/255,211/255],'filled')
            hold on
            
            for iL = 1:length(cumul_ID)
                if find(strcmpi(cumul_ID{iL},[num2str(iSess) '_' num2str(iAAC)]))
                    scatter(cumul_thetaphase_aac(iL),  cumul_ripmod_aac(iL),dotSize,'filled','m')
                end
            end
            xlabel('theta phase')
            xlim([-pi pi])
            ylabel('ripple mod')
            yline(1,':')
            
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
            LatencyFirstSpike = zeros(1,length(pulseEpochs));
            for iPulse = 1:length(pulseEpochs)
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
            
            ph_bin = linspace(-pi,pi,16);
            k = gaussian2Dfilter([10 10],[.5 .5]);
            
            
            imagesc(ph_mod.ph_bin,[],nanconvn((ph_mod.ph_rate(:,1:end-1,iAAC)),k),[min(linearize(ph_mod.ph_rate(:,1:end-1,iAAC)))...
                max(linearize(ph_mod.ph_rate(:,1:end-1,iAAC)))])
            
            hold on
            colormap('jet')
            plot(ph_mod.ph_bin,10+cos(ph_mod.ph_bin)*10,'w')
            set(gca,'ytick',0:10:length(ph_mod.freq),'yticklabel',round(ph_mod.freq(1:10:end)))
            set(gca,'ydir','normal')
            ylabel('Frequency (logscale)')
            xlabel('Theta Phase')
        end
        title('Phase Portrait')
        
        %%
        % % % % % % % % % % % % %
        % % TBD
        % % % % % % % % % % % % %
        subplot(6,4,16)
        
        
        %%
        % % % % % % % % % % % % %
        % %Spikes per Ripple Dist
        % % % % % % % % % % % % %
        
        subplot(6,4,21)
        if exist([basename '.ripspikes.analysis.mat'])
            load([basename '.ripspikes.analysis.mat'])
        else
            
            [spikesRipNum, numSpkPerCycPerRip] = getNumSpkRip(basepath,'units',aacs,'saveMat',true);
        end
        
        maxHistoRip = max(spikesRipNum{iAAC});
        edgesRip = 0:maxHistoRip;
        xtRip = [0.5:5:maxHistoRip+0.5];
        histogram(spikesRipNum{iAAC},edgesRip)
        xlabel('Num Spk in Rip')
        ylabel('count')
        set(gca,'YScale','log','XTick',xtRip-0.5,'XTickLabel', num2cell(xtRip-0.5))
        box off
        
        %%
        % % % % % % % % % % % % %
        % %Ripple cycle spikes
        % % % % % % % % % % % % %
        subplot(6,4,22)
        
        
        edgesRip = 0:1:5;
        histogram(cell2mat(numSpkPerCycPerRip{iAAC}),edgesRip)
        xlabel('Avg Number of Spikes per Ripple Cycle')
        ylabel('count')
        set(gca,'YScale','log','XTick',edgesRip,'XTickLabel', num2cell(edgesRip))
        box off
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
        
        sgtitle({['Unit ' num2str(iSess) '_' num2str(spikes.UID(iAAC))]},'interpreter','none')
        
        if doSave
            unitStr = ['D:\Data\SummaryPlot_AAC_' num2str(iSess) '_' num2str(iAAC)];
            savefig(gcf,[unitStr '.fig'])
            print(gcf,[unitStr '.pdf'],'-dpdf','-bestfit')
            %         append_pdfs(['E:\Dropbox\PD_Hpc\Progress\AAC\AAC_SummaryPlots_andCCG\SummaryPlot_AACs_20201218_speed2cms.pdf'],[unitStr '.pdf'])
            append_pdfs(['E:\Dropbox\PD_Hpc\RerunSummaryPlots_AACs_20210302_speed2cms.pdf'],[unitStr '.pdf'])
            
            delete([unitStr '.pdf'])
            close gcf
        end
    end
end


