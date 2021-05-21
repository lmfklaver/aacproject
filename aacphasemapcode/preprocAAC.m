  %Preprocess
  basepath = cd;, basename = bz_BasenameFromBasepath(basepath)
  
  addpath(genpath('C:\Users\lklaver\Documents\GitHub\CellExplorer\'))
session = sessionTemplate(basepath,'showGUI',true);
cell_metrics = ProcessCellMetrics('session', session);
rmpath(genpath('C:\Users\lklaver\Documents\GitHub\CellExplorer\'))
            out = runZetaStats20_100()

  %%
  % session params
    sessionInfo = bz_getSessionInfo(cd);
    
    params.nChans       = sessionInfo.nChannels;
    params.sampFreq     = sessionInfo.rates.wideband;
    params.Probe0idx    = sessionInfo.channels;
    clear aacs spikesAACs

    spikes = bz_GetSpikes;
    [~, ~, aacs] = splitCellTypes(basepath);
    
  
        for iAAC = 1:length(aacs)
            spikesAACs.times{iAAC} = spikes.times{aacs(iAAC)};
            spikesAACs.rawWaveform{iAAC} = spikes.rawWaveform{aacs(iAAC)};
            spikesAACs.UID(iAAC) = spikes.UID(aacs(iAAC));
        end
        
        
        load([basename '.optoStim.manipulation.mat'])
        load([basename '.cell_metrics.cellinfo.mat'])
        
        % saving
        opts.doSave         = 1;
        opts.doSaveFig      = 1;
        opts.saveMat = true;
        opts.doPlot = true;
        
        % ccg
        opts.ccgBinSize  = 0.001;
        opts.ccgDur      = 0.2;
        
        %getCCG OUT pulse
        [status, interval] = cellfun(@(a) InIntervals(a,optoStim.timestamps),spikesAACs.times, 'uni',false);
        
        for iAAC = 1:length(aacs);
            spikesAACs.nopulse{iAAC} = spikesAACs.times{iAAC}(~status{iAAC});
            spikesAACs.CE_CellType{iAAC} = cell_metrics.putativeCellType{aacs(iAAC)};
            spikesAACs.CE_FR(iAAC) =  cell_metrics.firingRate(aacs(iAAC));
            spikesAACs.CE_optmod(iAAC) = cell_metrics.optoStim_modulationIndex(aacs(iAAC));
        end
        
        [ccg,t]=CCG(spikesAACs.nopulse,[],'Fs',params.sampFreq, 'binSize',opts.ccgBinSize,'duration', opts.ccgDur, 'norm', 'rate');
        save([basename '.ccgaac.mat'],'ccg','t')
        
        
            STP = ShortTermPlasticity_axax_Kaiser_edit(basepath,'excitation');
            getPhaseMap_Kaiser_edit(basepath)