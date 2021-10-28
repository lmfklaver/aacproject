% % % AAC Project % % % 

%%%%%%% Session IN %%%%%%


% % % % % % % % % % % % %
% % Kilosort
% % % % % % % % % % % % %
addpath(genpath('C:\Users\lklaver\Documents\GitHub\Kilosort2'))
addpath(genpath('C:\Users\lklaver\Documents\GitHub\npy-matlab'))

cd 'C:\Users\lklaver\Documents\GitHub\Kilosort2'
% select preprocess master_ file and adapt config accordingly

%%

% % % % % % % % % % % % %
% % Spikesort
% % % % % % % % % % % % %

% in phy2

%% 

% % % % % % % % % % % % %
% % Validate Spikesorting
% % % % % % % % % % % % %

% Make CCGs and validate spikesorting - merge where necessary



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

% analogin channels 0-based
params.analoginCh.pulse     = 3;
params.analoginCh.wheel     = 1;
params.analoginCh.reward    = 0;

% analogin chR
params.analoginCh.pulse     = 1;
params.analoginCh.wheel     = 0;
params.analoginCh.reward    = 2;

%mouse 4
params.analoginCh.pulse     = 3;
params.analoginCh.wheel     = 'none'; %none 
params.analoginCh.reward    = 'none'; % none

%mouse 5
params.analoginCh.pulse     = 2;
params.analoginCh.wheel     = 'none'; %none 
params.analoginCh.reward    = 'none'; % none

%mouse 6
params.analoginCh.pulse     = 3;
params.analoginCh.wheel     = 7; 
params.analoginCh.reward    = 'none';

%%
% Ripple channel 0 based

params.RippleChan = 22;  %% NB: MANUALLY SELECT RIP CHANNEL FROM DAT BECAUSE FINDBESTRIPCHAN SOMETIMES SUCKS
%%

% % % % % % % % % % % % %
% % Set Analyses options
% % % % % % % % % % % % %

% options 

%%
% % % % % % % % % % % % %
% % Detect Pulses
% % % % % % % % % % % % %

cd(analogin_path)

% % First get info of analogin
analogin_file   = [basename, '_analogin.dat'];

[analogin] = getAnaloginVals(basepath,'wheelChan',params.analoginCh.wheel  ,'pulseChan',params.analoginCh.pulse);


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

%%
% % % % % % % % % % % % %
% % Make optoStim.manipulation.mat and pulses.events.mat
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
% % Preprocess LFP
% % % % % % % % % % % % %
checkLFPfile = dir('*.lfp');

if isempty(checkLFPfile)
    bz_LFPfromDat;
end

lfp = bz_GetLFP('all');

%%
% % % % % % % % % % % % %
% % Verify Ripple Chan
% % % % % % % % % % % % %
if exist([basename '.ripples.events.mat'],'file')
    load('mouse1_180412_2.ripples.events.mat')
    rippleChan = ripples.detectorinfo.detectionchannel;
else 
    rippleChan =params.RippleChan;
end




%% 



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
