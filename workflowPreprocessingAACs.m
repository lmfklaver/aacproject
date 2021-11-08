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

% Convert to Klusters (Make Clu Res Fet Spk files)

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
params.analoginCh.pulse     = 4;
params.analoginCh.wheel     = 2;
params.analoginCh.reward    = 1;

params.RippleChan = 0;%% NB: MANUALLY SELECT RIP CHANNEL FROM DAT BECAUSE FINDBESTRIPCHAN SOMETIMES SUCKS


% % % % % % % % % % % % %
% % Set Analyses options
% % % % % % % % % % % % %

% options 



% % % % % % % % % % % % %
% % Preprocess LFP
% % % % % % % % % % % % %
checkLFPfile = dir('*.lfp');

if isempty(checkLFPfile)
    bz_LFPfromDat;
end

lfp = GetLFP('all');


% % % % % % % % % % % % %
% % Detect Ripples
% % % % % % % % % % % % %

chanRip = params.RippleChan; 
[ripples] = bz_FindRipples(cd,chanRip,'durations',[50 150],...
    'thresholds',[1 3], 'passband',[100 250], 'EMGThresh', 0.95,'saveMat',true);

% Make sure you indeed have the highest ripple channel 
edit findPyramidalLayer.m

% Redo bz_FindRipples if layer is not the same

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

if params.saveMat
    save([basename '_pulseEpochs'], 'pulseEpochs')
end

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
optoStim.detectorinfo       = 'getPulseEpochs';

save([basename '.optoStim.manipulation.mat'],'optoStim')

pulses = optoStim;
save([basename '.pulses.events.mat'],'pulses');

% % % % % % % % % % % % %
% % Load Spikes
% % % % % % % % % % % % %

%%TOGGLE%%
% spikes = bz_LoadPhy; % if from phy output
spikes = bz_GetSpikes; % if from CluResFetSpk


% % % % % % % % % % % % %
% % Run CellExplorer
% % % % % % % % % % % % %

session = sessionTemplate(basepath,'showGUI',true);
% Make sure that all of this session info is correct 

% Remove buzcode from your path
rmpath(genpath(tsnebuzPath));
cell_metrics = ProcessCellMetrics('session', session,'showGUI',true); 
% make sure it's all correct, and to also select "other metrics" to make
% sure the optostim.manipulation.mat is excluded from calculating
% burstiness etc. 


% % % % % % % % % % % % %
% % Run ripmod code
% % % % % % % % % % % % %


% % % % % % % % % % % % %
% % Run STP code
% % % % % % % % % % % % %


% % % % % % % % % % % % %
% % Run getPhaseMap code
% % % % % % % % % % % % %




