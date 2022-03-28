function [spkInstPhaseStruct]=spkInstPhaseOutPulse(basepath, icell)
% This function analyzes the phase locking of each cell to the various
% frequencies wanting to be analyzed
% Rate is higher when at certain phase in certain frequencies
%
%   USAGE
%   [binnedHisto,histoEdges, spkInstPhaseStruct]=spkInstPhaseOutPulse(basepath, run, icell)
%
%   Dependencies: Buzcode, Englishlab\utilities
%
%   INPUTS
%   Name-value paired inputs:
%   'basepath'      - folder in which .STP.mat and 'phmod.mat' can be found
%                       (Required input, Default is pwd)
%   'basename'      - basefile name to load
%   'nfreq'         - How many frequencies to look at for phasemap
%                       (Default: 25)
%   'freqRange'     - (Default: [1 250])
%   'freqspace'     - 'log','lin';
%   'saveMat'       - save output yes or no
%   'muaThres'      - MUA phase-locking threshold to verify that there is
%                   enough power (Default: 0.05)
%   'nPhaseBins'    - (Default:32)
%   'epochs'        - allows you to only take spikes in (e.g. run)
%                       epochs
%
%   OUTPUTS
%   histogram
%   ph_portrait     - struct containing the following fields:
%   .ph_rate        - freqbin x phasebins x unit of spikerate
%   .ph_map         - struct with trans_phase and time
%   .ph_bin         - phasebins
%   .freq           -  frequencies over which phasemap is calculated
%   .nfreq          -  number of frequency bins
%
%   EXAMPLES
%  [binnedHisto,histoEdges, spkInstPhaseStruct]=spkInstPhaseOutPulse(basepath, cell)
%
%
%   HISTORY
%   Original code taken from Lianne Klaver
%   Edited by Earl Gilbert 2/2022



muaThres            = 0.05;


cd(basepath);

%% load output from short term plasticity analysis

basename = bz_BasenameFromBasepath(basepath);

load([basename  '.STP.mat']);
load([basename '.optoStim.manipulation.mat']);
load([basename '.ripples.events.mat'],'ripples');
load([basename '.run.states.mat']);
ripChan      = ripples.detectorinfo.detectionchannel;
[gd_eps]     =  get_gd_eps(basepath);


% pulls the channel from the ripples and loads the xml file
xmln = [basepath '/' basename '.xml'];
fname = [basepath '/' basename '.lfp'];
xml = LoadXml(xmln);

% Load in the LFP during good epochs (gd_eps) outside of stim periods
d   =[];
ts  =[];

for k = 1:size(gd_eps,1);
    dt = LoadBinary(fname,'nchannels',xml.nChannels,'channels',ripChan+1,'frequency',xml.lfpSampleRate,'start',gd_eps(k,1),'duration',diff(gd_eps(k,:)));
    d = [d;dt]; %concat all gd_eps
    ts = [ts ; [0:length(dt)-1]'/xml.lfpSampleRate + gd_eps(k,1)];
end


%load spikes
load([basename  '.spikes.cellinfo.mat']);

%define the number of frequencies with which to look at phase
%preference
freqband = [5 10]; % [5 10] % average phase is calculated over the entire [5 10]
% nfreq   = 25; % was echt 200 of zo LK % set number of frequencies
% freq    = logspace(log10(1),log10(250),nfreq); % set frequencies along a log scale
ph_bin  = linspace(-pi,pi,32); % set phases into 32 bins
% ph_rate = nan(nfreq-1,32,length(spkInstPhaseStruct.times)); % reserve ph_rate space for the number of frequencies and how they are binned


%get spike times OUT pulse
[status, interval] = cellfun(@(a) InIntervals(a,gd_eps),spikes.times, 'uni',false);           
for iUnit = 1:length(spikes.times);
    spikes.nopulse{iUnit} = spikes.times{iUnit}(status{iUnit});
end 

%make run.epochs an input
[status, interval] = cellfun(@(a) InIntervals(a, run.epochs),spikes.nopulse, 'uni',false);           
for iUnit = icell;
    spikes.runandgood{iUnit} = spikes.times{iUnit}(status{iUnit});
end

spkInstPhaseStruct.times = {spikes.runandgood{icell}};
spkInstPhaseStruct.times{1}     = sort(cell2mat(spkInstPhaseStruct.times')); % overbodig, want maar 1 cell - puts all cells from the recording into one variable I am not sure if this is a good idea
spkInstPhaseStruct.UID{1}  = 1;
%%
% %It's important for the multishank probes for sure
% Mabe check size(d) , how many chans. Maybe indeed keep the AAC in,
% rewrite the iUnit = 1:length, to iUnit = iAAC. maxWaveformCh of iiAAC ,
% and only load that in and use as (d)
%sweet
%cross check with Bandpassfilter in getPhasePortrait / getPhasePref
allfreq=freqband(1):freqband(2);
for b=1:(length(allfreq)-1);
    selfreqband=[allfreq(b) allfreq(b+1)];
    filtered=[];
%filter LFP
filtered.data = BandpassFilter(double(d), 1250 , selfreqband); % Filter the data for the specified frequency made above%get inst phaes and amplitude
filtered.phase  = InstPhase(filtered.data); % pull the instantaneous phase using the hilbert transform
filtered.amp    = InstAmplitude(filtered.data); % pull the instantaneous amplitude using hilbert transorm
filtered.timestamps     = ts ; % set timestamps of the filtered signal
filtered.sampleRate   = xml.lfpSampleRate; % sampling rate
filtered.filterparms.passband   = selfreqband;
PowerPhaseRateStruct.times{1} = cell2mat(spkInstPhaseStruct.times(1));
PowerPhaseRateStruct.UID = cell2mat(spkInstPhaseStruct.UID(1));
[PowerPhaseRatemap,~] = bz_PowerPhaseRatemap_sam(PowerPhaseRateStruct,filtered,'powernorm','none'); 
%%See if thresholding fixes it?
kp_pwr_bns = [];
kp_pwr_bns = PowerPhaseRatemap.powerbins(PowerPhaseRatemap.powerskew>muaThres);
    
        % Below here needs to be fixed for looping over many recordings
        if any(kp_pwr_bns)
        kp_pwr_bns = kp_pwr_bns(1); %keeping the first value from the modulated power
        
        %convert inst amplitude to normalized amplitude (as was done in
        %bz_PowerPhaseRatemap)
        
        amp = NormToInt(log10(filtered.amp),'Z',[0 Inf],filtered.sampleRate);
        
        filtered.phase(amp<kp_pwr_bns) = nan;
        end
        

%%

%Get Power/Phase at each spike
spkInstPhaseStruct.amp(b,:) =cell2mat(cellfun(@(X) interp1(filtered.timestamps,filtered.amp,X,'nearest'),...
    spkInstPhaseStruct.times,'uniformoutput',false));
spkInstPhaseStruct.phase(b,:) = cell2mat(cellfun(@(X) interp1(filtered.timestamps,filtered.phase,X,'nearest'),...
    spkInstPhaseStruct.times,'uniformoutput',false));
spkInstPhaseStruct.ph_bin = ph_bin;
histoEdges = ph_bin;
spkInstPhaseStruct.freqband(b,:)=selfreqband;
spkInstPhaseStruct.binnedHisto(b,:) = histcounts(spkInstPhaseStruct.phase(b,:),histoEdges);
end
spkInstPhaseStruct.meanbinnedHisto=mean(spkInstPhaseStruct.binnedHisto);


