%%%%%%% Session IN %%%%%%


% % % % % % % % % % % % %
% % Kilosort
% % % % % % % % % % % % %
addpath(genpath('C:\Users\lklaver\Documents\GitHub\Kilosort2'))
addpath(genpath('C:\Users\lklaver\Documents\GitHub\npy-matlab'))

cd 'C:\Users\lklaver\Documents\GitHub\Kilosort2'
% select preprocess master_ file and adapt config accordingly



% % % % % % % % % % % % %
% % Spikesort
% % % % % % % % % % % % %

% in phy2

%%

% % % % % % % % % % % % %
% % Set Session Specific Parameters
% % % % % % % % % % % % %

sessionInfo = bz_getSessionInfo(cd);

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

% analogin chR
params.analoginCh.pulse     = 2;
params.analoginCh.wheel     = 1;
params.analoginCh.reward    = 3;

%mouse 4
params.analoginCh.pulse     = 4;
params.analoginCh.wheel     = 1; %none 
params.analoginCh.reward    = 1; % none

%mouse 5
params.analoginCh.pulse     = 3;
params.analoginCh.wheel     = 1; %none 
params.analoginCh.reward    = 1; % none

%mouse 6
params.analoginCh.pulse     = 4;
params.analoginCh.wheel     = 8; 
params.analoginCh.reward    = 1; % none


params.RippleChan = 22;%% NB: MANUALLY SELECT RIP CHANNEL FROM DAT BECAUSE FINDBESTRIPCHAN SOMETIMES SUCKS


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

lfp = bz_GetLFP('all');


% % % % % % % % % % % % %
% % Detect Ripples
% % % % % % % % % % % % %

chanRip = params.RippleChan; 
[ripples] = bz_FindRipples(cd,chanRip,'durations',[50 150],...
    'thresholds',[1 3], 'passband',[100 250], 'EMGThresh', 0.95,'saveMat',true);

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

rhdfilename = 'info.rhd';
read_Intan_RHD2000_file_noprompt(rhdfilename)

analogin_file   = [basename, '_analogin.dat'];
[analogin.pulse, analogin.pos, analogin.reward, analogin.ts] = getAnaloginVals(basename,params,board_adc_channels,params);

% if params.saveMat
    save([basename '_analogin'], 'analogin','-v7.3')
% end

[pulseEpochs] = getPulseTimes(analogin); %% NB pulsethreshold needs to be set per recording

% if params.saveMat
    save([basename '_pulseEpochs'], 'pulseEpochs')
% end

checkPulEvt = dir('*evt.ait*');

if isempty(checkPulEvt)
    makePulseFile
end


% % % % % % % % % % % % %
% % Make optoManipulation
% % % % % % % % % % % % %

% [pulses] = makePulsesStruct(basename,pulseEpochs)

pulses.timestamps  = pulseEpochs;
pulses.peaks            = (pulseEpochs(:,2)-pulseEpochs(:,1))/2;
pulses.amplitude        = [];
pulses.amplitudeUnits   = [];
pulses.eventID          = ones(length(pulseEpochs),1);
pulses.eventIDlabels    = cell(1,length(pulses.timestamps));
pulses.eventIDlabels(:) = {'OptoStim - Arch'};
% pulses.eventIDlabels: cell array with labels for classifying various event types defined in stimID (cell array, Px1).
% pulses.eventIDbinary: boolean specifying if eventID should be read as binary values (default: false).
pulses.center           = pulses.peaks;%;
pulses.duration         = pulseEpochs(:,2)-pulseEpochs(:,1);
pulses.detectorinfo     = 'getPulseEpochs'; 

save(strcat(basename, '.pulses.events.mat'), 'pulses')

optoStim.timestamps = pulses.timestamps;
optoStim.stimID     = pulses.eventID;%(==2)
optoStim.center     = pulses.center;
optoStim.duration   = pulses.duration;

save([basename '.optoStim.manipulation.mat'],'optoStim')



% % % % % % % % % % % % %
% % Load Spikes
% % % % % % % % % % % % %

%     spikes = bz_GetSpikes;
spikes = bz_LoadPhy;

% % % % % % % % % % % % %
% % run Stats
% % % % % % % % % % % % %
runZetaStats20_100

% % % % % % % % % % % % %
% % Run CellExplorer
% % % % % % % % % % % % %

session = sessionTemplate(basepath,'showGUI',true);
cell_metrics = ProcessCellMetrics('session', session);

% @ Kaiser, Please add in any additional options about phase that you have \
% - I don't have those implemented yet. 



% % % % % % % % % % % % %
% % Run optmod code
% % % % % % % % % % % % %

% Lianne will do this


% % % % % % % % % % % % %
% % Run STP code
% % % % % % % % % % % % %

% Lianne will do this

% % % % % % % % % % % % %
% % Run getPhaseMap code
% % % % % % % % % % % % %

% Lianne will do this

