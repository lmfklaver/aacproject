%%%%%%% Session IN %%%%%%

buzcodePath = 'C:\Users\lklaver\Documents\GitHub\buzcode\';
tsnebuzPath = 'C:\Users\lklaver\Documents\GitHub\buzcode\externalPackages\tSNE_matlab\';
cellExplorerPath = 'C:\Users\lklaver\Documents\GitHub\CellExplorer\';
aacprojectPath = 'C:\Users\lklaver\Documents\GitHub\aacproject\';


addpath(genpath(buzcodePath))
addpath(genpath(aacprojectPath))


% % % % % % % % % % % % %
% % Kilosort
% % % % % % % % % % % % %
addpath(genpath('C:\Users\lklaver\Documents\GitHub\Kilosort2'))
addpath(genpath('C:\Users\lklaver\Documents\GitHub\npy-matlab'))

cd 'C:\Users\lklaver\Documents\GitHub\Kilosort2'
% select preprocess master_ file and adapt config accordingly


% Run through KS1


% % % % % % % % % % % % %
% % Spikesort
% % % % % % % % % % % % %

basepath = cd; basename = bz_BasenameFromBasepath(basepath);
load('rez.mat')
%because these rez.ops paths are the wrong paths sometimes, make sure to
%change them before running convert

ks1path = 'F:\AACRecordingsSorted\m231_201120_094939\Kilosort_2021-12-14_075457';

rez.ops.basepath = basepath; 
rez.ops.basename = basename; 
rez.ops.savepath = fullfile(ks1path,'ClusResFetSpk');
% rez.ops.savepath = [ks1path filesep 'ClusResFetSpk'];

save('rez.mat', 'rez')

if ~exist(rez.ops.savepath)
    mkdir(rez.ops.savepath)
end



% Convert to Klusters (Make Clu Res Fet Spk files) from
% KilosortWrapper-master;

ConvertKilosort2Neurosuite_KSW(rez)

% Copy XML into folder  (Provided that Klusters overwrites the "Units"
% section of the XML
% Error - no spike groups - copy anatomical groups to spike groups, 48/24/3
% Alternatively: Copy CluResFet into basepath first and then continue from
% there. 

% Manually sort in Klusters (from basepath preferably). 

% Or in phy2)

%%

% % % % % % % % % % % % %
% % Set Session Specific Parameters
% % % % % % % % % % % % %

sessionInfo         = bz_getSessionInfo(cd);

params.nChans       = sessionInfo.nChannels;
params.sampFreq     = sessionInfo.rates.wideband;
params.Probe0idx    = sessionInfo.channels;

% wheel params
% % params.radiusDisk   = 26; % in cm
% % params.circDisk     = 2*pi*params.radiusDisk;

% analogin channels
params.analoginCh.pulse     = 3;
params.analoginCh.wheel     = 2;
params.analoginCh.reward    = 1;

params.RippleChan = 23;%% NB: MANUALLY SELECT RIP CHANNEL FROM DAT BECAUSE FINDBESTRIPCHAN SOMETIMES SUCKS


% % % % % % % % % % % % %
% % Set Analyses options
% % % % % % % % % % % % %

% options 
% toggle


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

chanRip = params.RippleChan; 
[ripples] = bz_FindRipples(cd,chanRip,'durations',[50 150],...
    'thresholds',[1 3], 'passband',[100 250], 'EMGThresh', 0.95,'saveMat',true);

% Make sure you indeed have the highest ripple channel 
% edit findPyramidalLayer.m
edit findRippleLayerChan

% Redo bz_FindRipples if layer is not the same 
% Note Well: If different ripple layer per shank - > multiple ripples.events.mats?

% Check if Ripples are detected well
checkRipEvtfile = dir('*evt.rip*');
if isempty(checkRipEvtfile)
    makeRipFile
end


% % % % % % % % % % % % %
% % Detect Pulses
% % % % % % % % % % % % %

cd(analogin_path)

% % First get info of analogin
rhdfilename = [basename '_info.rhd'];
read_Intan_RHD2000_file_noprompt(rhdfilename)

analogin_file   = [basename, '_analogin.dat'];
[analogin.pulse, analogin.pos, analogin.reward, analogin.ts] = getAnaloginVals(basename,params,board_adc_channels,params);

if params.saveMat
    save([basename '_analogin'], 'analogin')
end

[pulseEpochs] = getPulseTimes(analogin); %% NB pulsethreshold needs to be set per recording

% if params.saveMat
%     save([basename '_pulseEpochs'], 'pulseEpochs')
% end

checkPulEvt = dir('*evt.ait*');

if isempty(checkPulEvt)
    makePulseFile
end


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
% % Load Spikes
% % % % % % % % % % % % %

% First make sure your cluresfetspk are inside the basepath with the xml
% and datfile
% make sure you comment out the default of grabbing a ks1path in
% bz_GetSpikes


%%TOGGLE%%
% spikes = bz_LoadPhy; % if from phy output
spikes = bz_GetSpikes('sortingMethod','clu'); % if from CluResFetSpk 


% % % % % % % % % % % % %
% % Run CellExplorer
% % % % % % % % % % % % %

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
cell_metrics = ProcessCellMetrics('session', session,'showGUI',true); 
% make sure it's all correct, and to also select "other metrics" to make
% sure the optostim.manipulation.mat is excluded from calculating
% burstiness etc. 

% % % % % % % % % % % % %
% % gd_eps
% % % % % % % % % % % % %

gd_eps=get_gd_eps(basepath);

% % % % % % % % % % % % %
% % Cluster Quality
% % % % % % % % % % % % %
clusters = getClusterQuality(basepath)

% % % % % % % % % % % % %
% % CCG in out
% % % % % % % % % % % % %

%Ask Lianne
[ccginout] = getCCGinout(basepath, spikes, pulseEpochs) %gd_eps?

% % % % % % % % % % % % %
% % RUN thngs
% % % % % % % % % % % % %

%getVelocity
%getRunEpochs ( 5cm/s)
selRunEpochsIdx = run.epochs(:,2)-run.epochs(:,1) >3; % longer than three seconds
getRunEpochs
[vel] = getVelocity(analogin,'doFigure',false,'downsampleFactor',3000);
[run] = getRunEpochs(basepath,vel,'minRunSpeed',minRunSpeed,'saveMat',true,'saveAs','.run2cm.states.mat');


% % % % % % % % % % % % %
% % Pulse PETH
% % % % % % % % % % % % %
[peth] = getPETH_epochs(basepath,'epochs',pulseEpochs,'timwin',[-0.4 0.4], ...
               'binSize', 0.01, 'saveAs', '.pethPulse.mat')

% % % % % % % % % % % % %
% % RipCCG
% % % % % % % % % % % % %
[ripple_ccg_mac] = getRipCCG(basepath,spikes,'epochs',gd_eps,'ccgbin', 0.02,'ccgdur', .8);
[ripple_ccg_mic] = getRipCCG(basepath,spikes,'epochs',gd_eps,'ccgbin', 0.001,'ccgdur', .1);
[ripple_ccg] = getRipCCG(basepath,spikes,'epochs',gd_eps,'save_mat',true)

% % % % % % % % % % % % %
% % Run ripmod code
% % % % % % % % % % % % %

ripmod = getRipMod(basepath, spikes, 'epochs', gd_eps, 'ccg', ripple_ccg,'baseTime',[-0.4 -0.3],'baselineAroundPeak',[-.05 .05],'savemat',false)
%get stats

% % % % % % % % % % % % %
% % Get PHAMP
% % % % % % % % % % % % %
phaseFrequencyAmplitude %script needs to turn into function

% % % % % % % % % % % % %
% % Burstiness
% % % % % % % % % % % % %
burstinesstest  %script needs to turn into function

% % % % % % % % % % % % %
% % Spikes Per Ripple Cycle
% % % % % % % % % % % % %
spikesRipAndCycle %script needs to turn into function

% % % % % % % % % % % % %
% % Run STP code
% % % % % % % % % % % % %
ShortTermPlasticity_axax_bare %Lianne

% % % % % % % % % % % % %
% % Run getPhaseMap code
% % % % % % % % % % % % %



% % % % % % % % % % % % %
% % Summaryplots code with generating .mat files
% % % % % % % % % % % % %
summaryplotsAACs

