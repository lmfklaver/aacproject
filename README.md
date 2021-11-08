# aacproject

# Start
summaryplotsAACs are a good main script to look at the different functions being called for\

# Parameters chosen
Run is defined as >5 cm/s\
Run onset is defined as >2 cm/s

# Toolboxes used
Buzcode: <a href="https://github.com/buzsakilab/buzcode">https://github.com/buzsakilab/buzcode</a>\
CellExplorer:<a href="https://github.com/petersenpeter/CellExplorer">https://github.com/petersenpeter/CellExplorer</a>\
ZETA: <a href="https://www.biorxiv.org/content/10.1101/2020.09.24.311118v1">https://www.biorxiv.org/content/10.1101/2020.09.24.311118v1</a>\
<a href="https://github.com/JorritMontijn/ZETA">https://github.com/JorritMontijn/ZETA</a>\
Cite as: eLife 2021;10:e71969 DOI: 10.7554/eLife.71969

# .Mat files that are being loaded in

## basic
basename.ripples.events.mat through buzcode\
basename.optoStim.manipulation.mat through getPulseEpochs, through a wrapper pulseEpochsWrapperCellExplorer.m\
basename.spikes.cellinfo.mat through buzcode (bz_GetSpikes.m or bz_LoadPhy.m)\

## from CellExplorer
basename.mono_res.cellinfo.mat through CellExplorer\
basename.cell_metrics.cellinfo.mat through CellExplorer\


To obtain these:\
>> rmpath(genpath('~\Documents\GitHub\buzcode\externalPackages\tSNE_matlab')) % tSNE in buzcode clashes with tSNE of CE. \
>> cd datafolder\
>> session = sessionTemplate(basepath,'showGUI',true);\
>> cell_metrics = ProcessCellMetrics % GUI will open, select other_metrics to exclude optostim from burstiness index) \
% will make changes to spike struct to include a cluID


## stats
basename.pethzeta.stats.mat using runZeta.m and ZETA toolbox \

## phasemaps and STP
basename.STP.mat through ShortTermPlasticity_axax.m\
basename.ph_mod.mat through getPhaseMap.m\

## CCG analyses
basename.ccginout.analysis.mat through getCCGinout.m

## peths
basename.pulsepeth.analysis.mat through getPETH_epochs

## For Running:
basename.analogin.mat through getAnaloginVals\
basename.run2cm.states.mat through getRunEpochs.m (first run getVelocity.m)

