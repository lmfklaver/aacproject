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
basename = bz_BasenameFromBasepath(cd)

sessionInfo = bz_getSessionInfo(cd);

params.nChans       = sessionInfo.nChannels;
params.sampFreq     = sessionInfo.rates.wideband;
params.Probe0idx    = sessionInfo.channels;

% wheel params
% % params.radiusDisk   = 26; % in cm
% % params.circDisk     = 2*pi*params.radiusDisk;

% analogin channels base 1
params.analoginCh.pulse     = 4;
params.analoginCh.wheel     = 3;
params.analoginCh.reward    = 1;

params.RippleChan = 6;%% NB: MANUALLY SELECT RIP CHANNEL FROM DAT BECAUSE FINDBESTRIPCHAN SOMETIMES SUCKS

save([basename '_params.mat'], 'params')

% % % % % % % % % % % % %
% % Set Analyses options
% % % % % % % % % % % % %

% options 


%%
% % % % % % % % % % % % %
% % Preprocess LFP
% % % % % % % % % % % % %
checkLFPfile = dir('*.lfp');

if isempty(checkLFPfile)
    bz_LFPfromDat;
end

% lfp = bz_GetLFP('all');


% % % % % % % % % % % % %
% % Detect Ripples
% % % % % % % % % % % % %

chanRip = params.RippleChan; 
[ripples] = bz_FindRipples(cd,chanRip,'durations',[50 150],...
    'thresholds',[1 2], 'passband',[100 250], 'EMGThresh', 0.95,'saveMat',true);

checkRipEvtfile = dir('*evt.rip*');
if isempty(checkRipEvtfile)
    makeRipFile
end

%%
% % % % % % % % % % % % %
% % Detect Pulses
% % % % % % % % % % % % %

%cd(analogin_path)

% % First get info of analogin
% rhdfilename = [basename '_info.rhd'];
% read_Intan_RHD2000_file_noprompt(rhdfilename)

analogin_file   = [basename, '_analogin.dat'];
% [analogin.pulse, analogin.pos, analogin.reward, analogin.ts] = getAnaloginVals(basename,params,board_adc_channels,params);
v = double(bz_LoadBinary([basename '_analogin.dat'],'nChannels', 8, 'channels', params.analoginCh.pulse, 'precision', 'uint16'));
analogin.pulse = v * 0.000050354;
clear v
sr = 30000
analogin.ts = (1:length(analogin.pulse))/sr;
[pulseEpochs] = getPulseTimes(analogin);
%% make puls evt file
checkPulEvt = dir('*evt.ait*');
if isempty(checkPulEvt)
    makePulseFile
end
save([basename '_pulseEpochs'], 'pulseEpochs')
%%
% % % % % % % % % % % % %
% % load/save pos/reward
% % % % % % % % % % % % %

analogin = [];
analogin.pos = double(bz_LoadBinary([basename '_analogin.dat'],'nChannels', 8, 'channels', params.analoginCh.wheel, 'precision', 'uint16')) * 0.000050354;
analogin.reward = double(bz_LoadBinary([basename '_analogin.dat'],'nChannels', 8, 'channels', params.analoginCh.reward, 'precision', 'uint16')) * 0.000050354;

save([basename '_analogin'], 'analogin')


%%

% if params.saveMat
%     save([basename '_analogin'], 'analogin')
% end
% added below line to getAnaloginVals
% [pulseEpochs] = getPulseTimes(analogin); %% NB pulsethreshold needs to be set per recording


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


%%
% % % % % % % % % % % % %
% % Load Spikes
% % % % % % % % % % % % %

%     spikes = bz_GetSpikes;
spikes = bz_LoadPhy




%%
% % % % % % % % % % % % %
% % Run CellExplorer
% % % % % % % % % % % % %

session = sessionTemplate(cd,'showGUI',true);
% Crtl+I - forces use of the .xml probe layout


cell_metrics = ProcessCellMetrics('session', session);

cell_metrics = CellExplorer('metrics',cell_metrics);





%% attempt state scoring

EMGFromLFP = bz_EMGFromLFP(cd)
SleepState = SleepScoreMaster(cd,'ignoretime', [pulseEpochs(1,1) pulseEpochs(end,2)],'overwrite', true)% exclude the times of stimulation

%% make input struct for the state editor
lfp = bz_GetLFP([12 4 26]) % top channel, best rip chan, bottom channel only do 3 channels
x = ones(1,3);
inputData.rawEeg = mat2cell(double(lfp.data),length(lfp.data),x)
inputData.eegFS = 1
inputData.Chs = lfp.channels
inputData.MotionType = 'File'
inputData.motion = double(bz_LoadBinary([basename '_analogin.dat'],'nChannels', 8, 'channels', params.analoginCh.wheel, 'precision', 'uint16', 'downsample', 30000)) * 0.000050354; %state editor motion needs data in a one hz format
% inputData.motion(end) = []; % this may be needed if you run into an error
% on line 983
clear lfp


TheStateEditor(basename, inputData) % this is the manual state editor, use to check the automation of the sleepscoremaster
% for troublshooting with Lianne
% inputData.motion = downsample(analogin.pos,30000)
inputData.motion(end) = []; % delete basename.eegstates.mat before rerunning the editor
save([basename '_inputData.mat'], 'inputData')



% pos = double(bz_LoadBinary([basename '_analogin.dat'],'nChannels', 8, 'channels', params.analoginCh.wheel, 'precision', 'uint16', 'downsample', 24)) * 0.000050354;

% @ Kaiser, Please add in any additional options about phase that you have \
% - I don't have those implemented yet.

% % % % % % % % % % % % %
% % runEpochs
% % % % % % % % % % % % %
[selRunEpochs]=doPETHRun


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

