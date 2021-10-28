% Unc5B R01 Grant Figures

% Set paths
% Paths
addpath(genpath('E:\Dropbox\Code\Englishlab\AAC\'))
% addpath(genpath('E:\Dropbox\Code\toolbox\buzcode\'))
% addpath('E:\Dropbox\Code\intan')
% addpath('E:\Dropbox\Code\buzsupport')

% Set environment INH/EXC


% Datadir
basepath = 'D:\Data\Axoaxonic_Data_Lianne';

% Session
% basename = 'u19_200310_135409';
% basename = 'u19_200313_120452' ;
% basename = 'u19_200313_155505';     % nice one for initial UNC5B grant figures
%basename = sessions{iSess}
% basename = 'u21_200305_153604';
basename = 'm175_200821_151859_2';

% Session specific paths
analogin_path   = fullfile([basepath filesep basename]);
spike_path      = fullfile([basepath filesep basename]);

cd(analogin_path)
%% Parameters / Options
% session params
sessionInfo = bz_getSessionInfo(cd);

params.nChans       = sessionInfo.nChannels;
params.sampFreq     = sessionInfo.rates.wideband;
params.Probe0idx    = sessionInfo.channels;

% wheel params
% % params.radiusDisk   = 26; % in cm
% % params.circDisk     = 2*pi*params.radiusDisk;

% analogin
params.analoginCh.pulse     = 4;
params.analoginCh.wheel     = 2;
params.analoginCh.reward    = 1;



%  wheel stuff
% % opts.fastTrialsOnly = 1;

% saving
opts.doSave         = 1;
opts.doSaveFig      = 1;
opts.saveMat = true;
 opts.doPlot = 0;

% ccg
opts.ccgBinSize  = 0.001;
opts.ccgDur      = 0.2;

% ripples
if strcmpi(basename,'u19_200310_135409')
    opts.rippleChan = 9;
elseif strcmpi(basename,'u19_200313_155505')
opts.rippleChan     = 57;
elseif strcmpi(basename,'u19_200313_120452')
    opts.rippleChan = 38;
elseif strcmpi(basename,'u21_200305_153604')
    opts.rippleChan = 5;
elseif strcmpi(basename,'m175_200821_151859_2')
    opts.rippleChan = 0;
end

opts.minRipInt      = 50;
opts.maxRipDur      = 150;

%% Load analogin parameters.
cd(analogin_path)

% First get info of analogin
rhdfilename = [basename '_info.rhd'];
read_Intan_RHD2000_file_noprompt(rhdfilename)

analogin_file   = [basename, '_analogin.dat'];
[analogin.pulse, analogin.pos, analogin.reward, analogin.ts] = getAnaloginVals(basename,params,board_adc_channels,opts);

if opts.saveMat
save([basename '_analogin'], 'analogin')
end

%% Get Pulse Epochs
[pulseEpochs] = getPulseTimes(analogin);

if opts.saveMat
save([basename '_pulseEpochs'], 'pulseEpochs')
end

checkPulEvt = dir('*evt.ait*');

if isempty(checkPulEvt)
   makePulseFile
end


%% Get Ripples and Ripple Stats

% check if lfp file present
checkLFPfile = dir('*.lfp');
if isempty(checkLFPfile)
    bz_LFPfromDat;
end

% before finding ripples remove the artifacts of pulse from datfile 
% RemoveArtefact_dat

% check if ripfile saved
checkRipEvt = dir('*ripples.events.mat*');
if ~isempty(checkRipEvt)
    load(checkRipEvt.name)
    chanRip = opts.rippleChan;
else
    lfp = bz_GetLFP('all');
    chanRip = opts.rippleChan;
    minRipInt = opts.minRipInt;
    maxRipDur = opts.maxRipDur;
    [ripples] = bz_FindRipples(cd,chanRip,'durations',[minRipInt maxRipDur],'thresholds',[1 3], 'EMGThresh', 0.95,'saveMat',opts.saveMat);
end

% to make an evt.pulse file for Neuroscope: makeRipFile
checkRipEvtfile = dir('*evt.rip*');
if isempty(checkRipEvtfile)
   makeRipFile
end
% 
% 
% % Get the RipStats
% lfp2Filt = bz_GetLFP(chanRip);
% lfpFiltBP = bandpass(double(lfp2Filt.data),[100 200],1250);
% timestamps = lfp2Filt.timestamps;
% [maps,data,stats] = bz_RippleStats(lfpFiltBP,timestamps,ripples,'durations',[-0.2 0.2]);
% 
% % RipplePeak in Pulse?
% [peakInPulse, pulseWithRip, aRP] = inh_getPeakInPulse(ripples.peaks, pulseEpochs);

%% Load the spikes

cd(spike_path)
spikes = bz_LoadPhy; %% Load in the clusters that are spike sorted with phy and are labeled 'good'

%% Obtain ccg and 't' for CCG

%getCCG
[ccg,t]=CCG(spikes.times,[],'Fs',params.sampFreq, 'binSize',opts.ccgBinSize,'duration', opts.ccgDur, 'norm', 'rate');

%plotCCG all units with 'k' reference unit
inh_CCG

%% New CCG
% ACG in and out STIM

%% Rasters centered on Pulse(Start)

timwin = [-1 1];
options.binSize = 0.01;
options.doPlot = true;
options.doSaveFig = false;

% all pulses
% [rate, count, time] = inh_rastersToPulse_new(spikes, pulseEpochs,timwin,ccg,t,opts, params);
[rate, count, time] = inh_rastersToPulse_new(spikes, pulseEpochs,timwin,ccg,t,options);

save rasterdata rate count time pulseEpochs timwin ccg t options spikes
%% Rasters centered on Ripple Peak

% all ripples - 4 panel figures 
% inh_rastersToRipple(spikes, ripplePeak,timwin,ccg,t)

% ripples inside (P) and outside (NP) pulse - 8 panel figures
% [rateP, countP, timeP, rateNP, countNP, timeNP] = inh_rastersToRipPeak_new(spikes, ripples.peaks,timwin,ccg,t, opts, params);
[rateP, countP, timeP, rateNP, countNP, timeNP] = inh_rastersToRipPeak_new(spikes, ripples.peaks,timwin,ccg,t, opts, params);
save rippledata rateP countP timeP rateNP countNP timeNP ripples timwin ccg t options spikes

%% Heatmaps centered on Pulse
% manualsort
% INTIdx =  [1 11 14 41 46];
% AACIdx = [3 13 61];
INTIdx = [];
AACIdx = [];
windowLength = 8; % pick sigma value for the gaussian, different per binSize -> figure this out

inh_heatmap_rates(rate, windowLength,INTIdx, AACIdx)

%% Gain Histograms
inh_calcGainRipple

% inh_calcGainPulse

%% Heatmaps centered on Ripple Peak
inh_heatmap_rates_ripples
inh_heatmap_rates_ripples_pulsenopulse

%% Ripple Stats Boxplots / Histograms
inh_rippleStats_pulsenopulse
