# aacproject

# Start
summaryplotsAACs are a good main script to look at the different functions being called for

# Parameters chosen
Run is defined as >5 cm/s
Run onset is defined as >2 cm/s

# .Mat files that are being loaded in

## basic
basename.ripples.events.mat through buzcode
basename.optoStim.manipulation.mat through getPulseEpochs, through a wrapper pulseEpochsWrapperCellExplorer.m
basename.spikes.cellinfo.mat through buzcode (bz_GetSpikes.m or bz_LoadPhy.m)
basename.mono_res.cellinfo.mat through CellExplorer
basename.cell_metrics.cellinfo.mat through CellExplorer

## stats
basename.pethzeta.stats.mat using FUNCTIONNAME.m and ZETA toolbox 


## phasemaps and STP
basename.STP.mat through ShortTermPlasticity_axax.m
basename.ph_mod.mat through getPhaseMap.m

## CCG analyses
basename.ccginout.analysis.mat through getCCGinout.m

## peths 
basename.pulsepeth.analysis.mat through getPETH_epochs

## For Running:
basename.analogin.mat through getAnaloginVals
basename.run2cm.states.mat through getRunEpochs.m (first run getVelocity.m)

