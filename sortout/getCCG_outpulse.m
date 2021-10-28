
       params.sampFreq = 30000; 
        % saving
        opts.doSave         = 1;
        opts.doSaveFig      = 1;
        opts.saveMat = true;
        opts.doPlot = true;
        
        % ccg
        opts.ccgBinSize  = 0.001;
        opts.ccgDur      = 0.2;
          
%getCCG OUT pulse
        [status, interval] = cellfun(@(a) InIntervals(a,pulseEpochs),spikes.times, 'uni',false);
        
        for iUnit= 1:length(spikes.times);
            spikes.nopulsetimes{iUnit} = spikes.times{iUnit}(~status{iUnit});
        end
        
        [ccg,t]=CCG(spikes.nopulsetimes,[],'Fs',params.sampFreq, 'binSize',opts.ccgBinSize,'duration', opts.ccgDur, 'norm', 'rate');
        
        
        %% plot ccg figure per unit
% if plotCCGperUnit
    for iUnit = 1:size(ccg,2)
        figure
        for nPlot = 1:size(ccg,3)
            subplot(7,6,nPlot)
            if iUnit == nPlot
            bar(t,ccg(:,iUnit,nPlot),'k')
            else
                bar(t,ccg(:,iUnit,nPlot))
            end
        end
    end
% end
