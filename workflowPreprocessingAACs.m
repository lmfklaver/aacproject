%%%%%%% Session IN %%%%%%

%This step is needed because there is a redundancy in the tSNE function in
%buzcode and in cell explorer, if you run this as a single script, you need
%to remove the buzcode tSNE function from path to successfully run cell
%explorer, also it adds some necessary packages
buzcodePath = 'C:\Users\English-Admin\OneDrive - Virginia Tech\Documents\GitHub\buzcode';
tsnebuzPath = 'C:\Users\English-Admin\OneDrive - Virginia Tech\Documents\GitHub\buzcode\externalPackages\tSNE_matlab';
cellExplorerPath = 'C:\Users\English-Admin\OneDrive - Virginia Tech\Documents\GitHub\CellExplorer';
aacprojectPath = 'C:\Users\English-Admin\OneDrive - Virginia Tech\Documents\GitHub\lklaver\aacproject';


addpath(genpath(buzcodePath))
addpath(genpath(aacprojectPath))

%Set basepath and basename per buzcode formatting
basepath = cd; basename = bz_BasenameFromBasepath(basepath);

% % % % % % % % % % % % %
% % Kilosort
% % % % % % % % % % % % %


% Run through KS1
% See C:\Users\English-Admin\OneDrive - Virginia Tech\EnglishLab\Users\Earl\Matlab\SpikeSorting_SuperFolder\KS1_Run_BatchTemp.m


% % % % % % % % % % % % %
% % Set Session Specific Parameters
% % % % % % % % % % % % %
%May or may not be necessary
sessionInfo         = bz_getSessionInfo(cd);

params.nChans       = sessionInfo.nChannels;
params.sampFreq     = sessionInfo.rates.wideband;
params.Probe0idx    = sessionInfo.channels;

% wheel params
% % params.radiusDisk   = 26; % in cm
% % params.circDisk     = 2*pi*params.radiusDisk;

% analogin channels
params.analoginCh.pulse     = 6;
params.analoginCh.wheel     = 1;
params.analoginCh.reward    = 9;

params.RippleChan = 20;%% NB: MANUALLY SELECT RIP CHANNEL FROM DAT BECAUSE FINDBESTRIPCHAN SOMETIMES SUCKS


% % % % % % % % % % % % %
% % Preprocess LFP
% % % % % % % % % % % % %
checkLFPfile = dir('*.lfp');

if isempty(checkLFPfile)
    bz_LFPfromDat;
end

lfp = bz_GetLFP('all');


% % % % % % % % % % % % %
% % Detect Ripples
% % % % % % % % % % % % %
chanRip = 57; 
[ripples] = bz_FindRipples(cd,chanRip,'durations',[30 250],...
    'thresholds',[.5 1], 'passband',[100 250], 'EMGThresh', 0.99,'saveMat',true);
%RippleStats from buzcode
ripplelfp = bz_GetLFP(ripples.detectorinfo.detectionparms.channel);
rippleFilt = bz_Filter(ripplelfp, 'passband', [100 250]);
rippleFilt.data = double(rippleFilt.data)*0.195;% convert to microvolts
[ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(rippleFilt.data,rippleFilt.timestamps,ripples);

%Custom code for sharp wave detection - probably not well written, can
%optimize
SWChan=26
ripples=findSWAmp(basepath,ripples,SWChan);
[ripples.SW] = findSharpWaves('ripples',ripples,'rippleChannel',ripples.detectorinfo.detectionparms.channel,'SWChannel',SWChan)
save([basename '.ripples.events.mat'],'ripples');

% Make sure you indeed have the highest ripple channel 
% edit findPyramidalLayer.m
% edit findRippleLayerChan

% Redo bz_FindRipples if layer is not the same 
% Note Well: If different ripple layer per shank - > multiple ripples.events.mats?

% Check if Ripples are detected well in NeuroScope2


% % % % % % % % % % % % %
% % Load Spikes
% % % % % % % % % % % % %

spikes=bz_LoadPhy_CellExplorer %_CE found in utilites from kaiser so CE will work

% % % % % % % % % % % % %
% % Detect Pulses
% % % % % % % % % % % % %

% % First get info of analogin
rhdfilename = [basename '_info.rhd'];
read_Intan_RHD2000_file_noprompt(rhdfilename)

analogin_file   = [basename, '_analogin.dat'];
%Depending on how experiment was done, you can assign numbers to pulseChan,
%wheelChan, etc to create analogin.pulse analogin.wheel
[analogin] = getAnaloginVals(basepath,'pulseChan',1); %Base 1


%Right now we use the digitalin on the intan for pulses, it's more reliable
%so this should work for future recording sessions
fileinfo = dir([basename '_digitalin.dat']);
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen([basename '_digitalin.dat'], 'r');
digital_word = fread(fid, num_samples, 'uint16')';
fclose(fid);
analogin.pulse = (bitand(digital_word, 2^6)); %2^digital in channel number, ie if 2^6, pulse was going to digitalin channel 6

%Sanity check that you're pulling the right channel
unique(analogin.pulse)
plot(analogin.pulse)

%Soft fix for hardcoding ts and sr to analogin struct
sr=30000
analogin.ts      = (1:length(analogin.pulse))/sr;
analogin.sr      = sr;
save([basename '_analogin.mat'],'analogin')


%This detects the pulse epochs, there is some manual manipulation you may
%need to do, read through this function carefully
[pulseEpochs] = getPulseTimes(analogin); %% NB pulsethreshold needs to be set per recording


% % % % % % % % % % % % %
% % Make optoManipulation
% % % % % % % % % % % % %

optoStim.timestamps         = pulseEpochs;
optoStim.peaks              = (pulseEpochs(:,2)-pulseEpochs(:,1))/2;
optoStim.amplitude          = [];
optoStim.amplitudeUnits     = [];
optoStim.eventID            = ones(length(pulseEpochs),1);
optoStim.eventIDlabels      = cell(length(optoStim.timestamps),1);
optoStim.eventIDlabels(:)   = {'OptoStim - ChR'};
optoStim.center             = optoStim.peaks;%;
optoStim.duration           = pulseEpochs(:,2)-pulseEpochs(:,1);
optoStim.detectorinfo       = 'getPulseTimes';

save([basename '.optoStim.manipulation.mat'],'optoStim')

pulses = optoStim;
save([basename '.pulses.events.mat'],'pulses');
save([basename '.optoStim.events.mat'],'optoStim');

% % % % % % % % % % % % %
% % RUN things
% % % % % % % % % % % % %
%Pulls out running data from analogin.pos
minRunSpeed = 10
minRunLength = 3
[vel] = getVelocity(analogin,'circDisk',2*pi*26,'doFigure',true,'downsampleFactor',3000);
[run] = getRunEpochs(basepath,vel,'minRunSpeed',minRunSpeed,'saveMat',true,'saveAs','.run.states.mat');
yval=[]
yval(1:length(run.epochs(:,1)))=40
plot(vel.laps.pos_in_cm)
hold on
scatter(run.index(1,:),yval,'k*')
scatter(run.index(2,:),yval,'r*')

% % % % % % % % % % % % %
% % Run Behavior for Theta Cell Explorer
% % % % % % % % % % % % %

% Will need to define theta channel in CE and build this struct
% basename.animal.behavior.mat where animal.speed and animal.time and sr

% % % % % % % % % % % % %
% % Run CellExplorer
% % % % % % % % % % % % %
basepath = cd
addpath(genpath(cellExplorerPath))
session = sessionTemplate(basepath,'showGUI',true);


% Make sure that all of this session info is correct 
% and append any session info you want to add: 
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Additional animal metadata (https://cellexplorer.org/datastructure/data-structure/)
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    
%     session.animal.probeImplants.probe = '' %name of probe implanted
%     session.animal.probeImplants.brainRegion = '' %brain region
%     session.animal.probeImplants.ap = []; %Anterior-Posterior coordinate (mm)
%     session.animal.probeImplants.ml = []; %Medial-Lateral coordinate (mm)
%     session.animal.probeImplants.depth = []; %: Implantation depth (mm)
%     session.animal.probeImplants.ap_angle = []; %ap-angle of probe implantation (degrees)
%     session.animal.probeImplants.ml_angle = []; % ml angle of probe implantation (degrees)
%     session.animal.probeImplants.rotation = []; % rotation of probe (degrees)

% % % load('chanMap.mat')
% % z=sortrows([chanMap ycoords xcoords],1);
% % chanCoords.x = z(:,3);
% % chanCoords.y = z(:,2);
% % chanCoords.verticalSpacing = 20;
% % session.extracellular.chanCoords.x = z(:,3);
% % session.extracellular.chanCoords.y = z(:,2);
% % session.extracellular.chanCoords.verticalSpacing = 20;
% % save([basename '.session.mat'], 'session')
% % save([basename '.chanCoords.channelInfo.mat'], 'chanCoords')
% %     

%     session.animal.opticFiberImplants.opticFiber = ''; % optic fiber implanted
%     session.animal.opticFiberImplants.brainRegion =''; %brain region
%     session.animal.opticFiberImplants.ap = '';% : Anterior-Posterior coordinate (mm)
%     session.animal.opticFiberImplants.ml = '';% : Medial-Lateral coordinate (mm)
%     session.animal.opticFiberImplants.depth = '';%: Implantation depth (mm)
%     session.animal.opticFiberImplants.ap_angle = '';%: ap-angle of probe implantation (degrees)
%     session.animal.opticFiberImplants.ml_angle = '';%: ml angle of probe implantation (degrees)
%     session.animal.opticFiberImplants.notes = '';%: notes

%     session.animal.surgeries.date = ''; %
%     session.animal.surgeries.start_time = ''; 
%     session.animal.surgeries.end_time = '';
%     session.animal.surgeries.weight = '';
%     session.animal.surgeries.type_of_surgery = 'Acute';
%     session.animal.surgeries.room = '';
%     session.animal.surgeries.persons_involved = '';
%     session.animal.surgeries.anesthesia = 'Isoflurane';
%     session.animal.surgeries.analgesics == '';
%     session.animal.surgeries.antibiotics = '';
%     session.animal.surgeries.complications = '';
%     session.animal.surgeries.notes = '';

%     session.animal.virusInjections.virus = 'ChR';
%     session.animal.virusInjections.brainRegion ='';
%     session.animal.virusInjections.injection_schema = '';
%     session.animal.virusInjections.injection_volume = '';
%     session.animal.virusInjections.injection_rate = '';
%     session.animal.virusInjections.ap = [];
%     session.animal.virusInjections.ml = [];
%     session.animal.virusInjections.depth = [];
%     session.animal.virusInjections.ap_angle = [];
%     session.animal.virusInjections.ml_angle = [];
%       save([basename '.session.mat'],'session')

% session.channelTags.Theta.channels = 64;              % Theta channel
% session.channelTags.Ripple.channels = params.RippleChan;             % Ripple channel
% session.channelTags.RippleNoise.channels = 1;         % Ripple Noise reference channel
% session.channelTags.Cortical.electrodeGroups = 3;     % Cortical spike groups
% session.channelTags.Bad.channels = 1;                 % Bad channels
% session.channelTags.Bad.electrodeGroups = 1;          % Bad spike groups (broken shanks)

% session.epochs{1}.name = '';
% session.epochs{1}.startTime =  0;
% session.epochs{1}.stopTime = [];
% session.epochs{1}.behavioralParadigm = '';
% session.epochs{1}.builtMaze = '';
% session.epochs{1}.mazeType = '';
% session.epochs{1}.manipulations ='';


% Remove buzcode from your path
rmpath(genpath(tsnebuzPath));
% Run CellExplorer, 
%NEED TO INCLUDE VERTICAL SPACING IN EXTRACELLULAR OR WILL FAIL
cell_metrics = ProcessCellMetrics('session', session,'showGUI',true);


% make sure it's all correct, and to also select "other metrics" to make
% sure the optostim.manipulation.mat is excluded from calculating
% burstiness etc. 

% % % % % % % % % % % % %
% % gd_eps
% % % % % % % % % % % % %

% Selects Epochs of time >x Seconds away from stimulation times. This is
% important for optotagging experiments because we want to identify cells
% through optotagging, but we want to study them outside of manipulation
gd_eps=get_gd_eps(basepath);


% % % % % % % % % % % % %
% % Theta Modulation Index
% % % % % % % % % % % % 
%This script only selects spikes that occured during thetaEpochs
basename=bz_BasenameFromBasepath(cd)
load([basename '.spikes.cellinfo.mat'])
load([basename '.run.states.mat'])
load([basename '.gd_eps.mat'])
sr=30000
selSpikes.times=[];
gdSpikes.times=[];
for j = 1:length(spikes.times)
    [status] = InIntervals(spikes.times{j},gd_eps);
    gdSpikes.times{j}=spikes.times{j}(status);
end
for j = 1:length(gdSpikes.times);
    [status] = InIntervals(gdSpikes.times{j},thetaEpochs.intervals);
    selSpikes.times{j}=gdSpikes.times{j}(status);
    selSpikes.total(j)=sum(size(selSpikes.times{j},1));
end
selSpikes.numcells=size(selSpikes.times,2);
acg_metrics = calc_ACG_metrics(selSpikes,sr);
thetaModulationIndex=acg_metrics.thetaModulationIndex;
save([basename '.thetamodulationindex.mat'],'thetaModulationIndex')
addpath(genpath(buzcodePath))




% % % % % % % % % % % % %
% % Cluster Quality
% % % % % % % % % % % % %

%Is a measure of quality of spike sorting
[clusters] = getClusterQuality(basepath)
save([basename '.ClusterQuality.analysis.mat'], 'clusters') ;

% % % % % % % % % % % % %
% % CCG in out
% % % % % % % % % % % % %
[pulseEpochs] = optoStim.timestamps;
[ccginout] = getCCGinout(basepath, spikes, pulseEpochs); %gd_eps?


% % % % % % % % % % % % %
% % Pulse PETH
% % % % % % % % % % % % %
% Peri-event time histogram aligned to start of stimulation

% Leftover from looking at stims that happened during theta epochs, maybe
% useful later but not relevant
%[status]=InIntervals(optoStim.timestamps(:,1),thetaEpochs.intervals);

[pulsepeth] = getPETH_epochs(basepath,'epochs',optoStim.timestamps(:,1),'timwin',[-1 1], ...
               'binSize', 0.01,'long',true);
save([basename '.pulsepeth.analysis.mat'], 'pulsepeth') ;

% % % % % % % % % % % % %
% % Ripple PETH
% % % % % % % % % % % % %  
% Peri-event time histogram aligned to peak of SPW-Rs that are contained in
% gd_eps

[status]=InIntervals(ripples.peaks,gd_eps);
gd_ripplepeaks=ripples.peaks(status);
[ripplepeth] = getPETH_epochs(basepath,'epochs',gd_ripplepeaks,'timwin',[-0.5 0.5], ...
                     'binSize', 0.01,'long',true);
save([basename '.ripplepeth.analysis.mat'], 'ripplepeth') ;

% % % % % % % % % % % % %
% % Ripple PETH
% % % % % % % % % % % % %
% Peri-event time histogram aligned to peak of SPW-Rs that are contained in
% optoStim

[status]=InIntervals(ripples.peaks,optoStim.timestamps);
gd_ripplepeaks=ripples.peaks(status);
[ripplepeth] = getPETH_epochs(basepath,'epochs',gd_ripplepeaks,'timwin',[-0.5 0.5], ...
                     'binSize', 0.01,'long',true);
save([basename '.ripplepethOpto.analysis.mat'], 'ripplepeth') ;

% % % % % % % % % % % % %
% % RipCCG
% % % % % % % % % % % % %
% Similar to Ripple PETH, but instead is analysis of ripple ccg (think of a
% ccg in spike sorting, but the ripple peaks are considered a spike) This
% is the same output as a CCG function, but END is the CCGs to the ripples

%[ripple_ccg_mac] = getRipCCG(basepath,spikes,'epochs',gd_eps,'ccgbin', 0.02,'ccgdur', .8);
%[ripple_ccg_mic] = getRipCCG(basepath,spikes,'epochs',gd_eps,'ccgbin', 0.001,'ccgdur', .1);
[ripple_ccg] = getRipCCGFixed(basepath,spikes,'epochs',gd_eps,'ccgbin', 0.001,'ccgdur', 1,'saveMat',false);
%[ripple_ccgSTIM] = getRipCCGSTIM(basepath,spikes,'ccgbin', 0.01,'ccgdur', 1,'saveMat',true);
save([basename '.ripple_ccg1ms.analysis.mat'], 'ripple_ccg')

% % % % % % % % % % % % %
% % CCGs of spikes only within Gd_time, no rips
% % % % % % % % % % % % %

[Nonrips] = ExcludeIntervals(gd_eps,ripples.timestamps);
[NonripspikeCCG] = getCCGinout(basepath, spikes, Nonrips)
save([basename '.NonripspikeCCG.analysis.mat'], 'NonripspikeCCG')

% % % % % % % % % % % % %
% % CCGs of spikes only outside of rips and stim
% % % % % % % % % % % % %
[status,interval]=InIntervals(ripples.peaks(:,1),gd_eps); 
ripstart=ripples.timestamps(:,1);
ripend=ripples.timestamps(:,2);
gdrips=[];
gdrips(:,1) = ripstart(status);
gdrips(:,2) = ripend(status);
[Congdrips] = ConsolidateIntervals(gdrips);
[ripspikeCCG] = getCCGinout(basepath, spikes, Congdrips)
save([basename '.ripspikeCCG.analysis.mat'], 'ripspikeCCG')

% % % % % % % % % % % % %
% % Run ripmod code %Defunct, use cell explorer instead
% % % % % % % % % % % % %
[ripmod] = getRipMod(basepath, spikes, 'epochs', gd_eps, 'ccg', ripple_ccg,'baseTime',[-0.4 -0.3],'baselineAroundPeak',[-.05 .05],'saveMat',true);
%get stats

[status,interval]=InIntervals(ripples.peaks(:,1),gd_eps); %Detect ripples outside of stim
ripstart=ripples.timestamps(:,1);
ripend=ripples.timestamps(:,2);
gdrips=[];
gdrips(:,1) = ripstart(status)-.05;
gdrips(:,2) = ripend(status)+.05;
[Congdrips] = ConsolidateIntervals(gdrips)
[SpikeLFPCouplingGdEps]=bz_GenSpikeLFPCoupling(spikes,lfp,'frange',[120 250],'nfreqs',1,'spikeLim',1000000,'cellclass',allcelltypes,'int',Congdrips,'channel',unique([ripples.detectorinfo.detectionparms.channel spikes.maxWaveformCh(aacs)]))
save([basename '.SpikeLFPCouplingGdEps1F.mat'],'SpikeLFPCouplingGdEps')

% % % % % % % % % % % % %
% % Burstiness % I think we just use cell explorer
% % % % % % % % % % % % %
[burstIndex] = burstinessMizuseki_epochs(basepath,spikes,'epochs',gd_eps,'saveMat',true)

% % % % % % % % % % % % %
% % Zeta
% % % % % % % % % % % % %
[zeta] = runZeta(basepath,optoStim.timestamps(:,1),'saveMat',true);
% % [zetaSTIMRips] = runZeta(basepath,ripspikes.ONrips.timestamps(:,1),spikes,'saveMat',true);
[zeta] = runZeta(basepath,ripspikes.OFFrips.timestamps(:,1),spikes,'saveMat',true,'saveAs','.ripzeta.stats.mat');

% % % % % % % % % % % % %
% % Cell Types
% % % % % % % % % % % % %
[pyrs, ints, aac] = splitCellTypes(basepath); %Change to ignore Cell explorer and just look at PETH?

% Manually Check PulsePETHs
for i=1:size(pulsepeth.rate,1)        
    afigure = figure,;
    hold on;
    %load([basename '.pulsepeth.analysis.mat']);
            h2 = histogram('BinEdges',pulsepeth.timeEdges, ...
                'BinCounts',pulsepeth.rate(i,:));
            box off
            title(['PETH Stim' num2str(i)]);
            xlabel('time(s)');
            ylabel('spikes/s');
            h2.EdgeColor = 'none';
            h2.FaceColor = 'k';
            xlim(pulsepeth.timwin);
            line([0 0], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

% Plot a line at 0.3
line([0.3 0.3], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

% Add labels and title
    hold off
end
for i = 1:length(pulsepeth.trials);
    plotSpkOffset = 0;
    selTrialsPulse = pulsepeth.trials{i};
    figure,;
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
            ylim([ylimits(1) plotSpkOffset]);
            xlim(pulsepeth.timwin);
            line([0 0], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

% Plot a line at 0.3
line([0.56 0.56], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

end
%%Manually enter AACs
aacs=[6 17 18 37 45];
for selac=1:length(aacs)
pyrsnotaacs = pyrs~=aacs(selac);
pyrs=pyrs(pyrsnotaacs);
end
for selac=1:length(aacs)
intsnotaacs = ints~=aacs(selac);
ints=ints(intsnotaacs);
end
allcelltypes = cell(1,size(spikes.times,2));
allcelltypes(ints)={'int'};
allcelltypes(pyrs)={'pyr'};
allcelltypes(aacs)={'aac'};
save([basename '_celltypes'],'aacs', 'pyrs', 'ints', 'allcelltypes');
% % % % % % % % % % % % %
% % Spikes Per Ripple Cycle
% % % % % % % % % % % % %

[ripspikes] = getNumSpkRip(basepath,'units','all','saveMat',false)
[ripspikesNoStim] = getNumSpkRipNoStim(basepath,'units','all','saveMat',false)
[ripspiketime] = getRipSpkTime(basepath,gd_eps,'units','all','saveMat',true)
[ripspiketimeSTIM] = getRipSpkTimeSTIM(basepath,'units','all','saveMat',true)
[ripSTA] = getRipSTA(basepath,gd_eps,'units','all','saveMat',true)
[ripSTASTIM] = getRipSTASTIM(basepath,'units','all','saveMat',true)
[ripple_ccg_ON] = getRipCCGFixed(basepath,spikes,'epochs',ripspikes.ONrips.timestamps,'ccgbin', 0.001,'ccgdur', 1,'saveMat',false);
[ripple_ccg_OFF] = getRipCCGFixed(basepath,spikes,'epochs',ripspikes.OFFrips.timestamps,'ccgbin', 0.001,'ccgdur', 1,'saveMat',false);

save([basename '.ripple_ccg_STIM.mat'],'ripple_ccg_ON','ripple_ccg_OFF')

% % % % % % % % % % % % %
% % Spikes in Ripple Polar Plot
% % % % % % % % % % % % %

[RipSpikePhase]=getRipSpikePhase(cd)
[pval m] = circ_rtest(RipSpikePhase{i})
circ_plot(RipSpikePhase{i},'hist',[],20,true,true,'linewidth',2,'color','r')
title('HPC theta phase of HOR peaks')
legend(['P-value = ' num2str(pval)],'Mean theta phase')



% % % % % % % % % % % % %
% % Spikes Rates During Run
% % % % % % % % % % % % %
[runspikes] = getNumSpkRun(basepath,'units','all','saveMat',true)


% % % % % %
% % Run STP code % We don't use
% % % % % % % % % % % % %
[STP] = ShortTermPlasticity(basepath,'saveMat',true);
% % % % % % % % % % % % %
% % Run getPhaseMap code
% % % % % % % % % % % % %
runepochs=load([basename '.run.states.mat']);
pulseEpochs = optoStim.timestamps;
[pulseinrun] = findPulseInRun(runepochs, pulseEpochs);
[gd_run,indices] = SubtractIntervals(runepochs.run.epochs,pulseEpochs);
stimrun=pulseEpochs(pulseinrun,:);
[ph_mod_stimrun] = getPhasePref(basepath, 'epochs', stimrun,'freqRange',[5 10],'saveMat',false);
[ph_mod_gdrun] = getPhasePref(basepath, 'epochs', gd_run,'freqRange',[5 10],'saveMat',false);

[thetaEpochs] = detectThetaEpochs('bandpass',[4 10],'powerThreshold',1.2);
[thetanorip,indices] = SubtractIntervals(thetaEpochs.intervals,ripples.timestamps);
thetaEpochs.thetanoripintervals=thetanorip
save([basename '.thetaEpochs.states.mat'],'thetaEpochs')
[ph_mod] = getPhasePref(basepath, 'epochs', thetanorip,'freqRange',[5 10],'saveMat',true);
[ph_mod] = getPhasePref(basepath, 'epochs', thetanorip,'freqRange',[39 50],'saveMat',true);
[ph_portrait] = getPhasePortrait(basepath, 'epochs', thetanorip,'saveMat',true);
[ph_portrait_rip] = getPhasePortrait(basepath,'saveMat',true,'saveAs','.ph_portrait_rip.analysis.mat');

%%Hippo CookBook

computePhaseModulation('excludeIntervals',optoStim.timestamps)

% % % % % % % % % % % % %
% % Rank Order Spike Times
% % % % % % % % % % % % %
bz_getRipSpikes('basepath', basepath, 'saveMat', true);
[eventIDs]=InIntervals(ripples.peaks,optoStim.timestamps);
[rankStats] = RankOrder('eventIDs',double(eventIDs));
figure,histogram(rankStats.rankClusters(~logical(eventIDs)),'Normalization','Probability')
title('ID Cluster Stim Ripples')
figure,
histogram(rankStats.rankClusters(logical(eventIDs)),'Normalization','Probability')
title('ID Cluster Control Ripples')


% % % % % % % % % % % % %
% % Linear Mixed Model
% % % % % % % % % % % % %

predictors = {'cumul_FR_aac', 'cumul_abratio_aac', 'cumul_burstiness_aac', ...
              'cumul_thetam_aac', 'cumul_ratepowercorr_aac', ...
              'cumul_ripratepowcorrpremono', 'cumul_thetar_aac'};

results = [];

for i = 1:length(predictors)
    predictor = predictors{i};
    
    % Keep rows where all required columns are not NaN
    valid_rows = ~isnan(T.cumul_ripmodCE_aac) & ...
                 ~isnan(T.(predictor)) & ...
                 ~isnan(T.Mouse_ID) & ...
                 ~isnan(T.Session_ID);
             
    T_clean = T(valid_rows, :);
    
    % Only run the model if there's enough data
    if height(T_clean) > 10
        formula = sprintf('cumul_ripmodCE_aac ~ %s + (1|Mouse_ID) + (1|Session_ID)', predictor);
        lme = fitlme(T_clean, formula);
        coef = lme.Coefficients;
        pval = coef.pValue(2); % p-value of the predictor
        rsq = lme.Rsquared.Adjusted;
        results = [results; {predictor, pval, rsq, height(T_clean)}];
    else
        % Not enough data to fit model
        results = [results; {predictor, NaN, NaN, height(T_clean)}];
    end
end

results_table = cell2table(results, 'VariableNames', ...
    {'Predictor', 'pValue', 'AdjRsquared', 'NumObservations'});
disp(results_table)  

% % % % % % % % % % % % %
% % Summaryplots code with generating .mat files
% % % % % % % % % % % % %
summaryplotsAACs

