function [binnedHisto,histoEdges, spkInstPhaseStruct]=spkInstPhaseOutPulse(basepath, iAAC)

%% load output from short term plasticity analysis

basename = bz_BasenameFromBasepath(basepath);

load([basename  '.STP.mat'])
load([basename '.optoStim.manipulation.mat'])


%% get correct LFP channel (one from ripple)

% get ripple
fils  = getAllExtFiles(basepath,'rip',1);
rip = LoadEvents(fils{1});

% pulls the channel from the ripples and loads the xml file
ch = str2num(rip.description{1}(regexp(rip.description{1},'[0-9]')));
xmln = [basepath '/' basename '.xml'];
fname = [basepath '/' basename '.lfp'];
xml = LoadXml(xmln);

% load lfp 
d =[];ts =[];
dt = LoadBinary(fname,'nchannels',xml.nChannels,'channels',ch+1,'frequency',xml.lfpSampleRate);
d  = [d;dt];
ts = [ts ; [0:length(dt)-1]'/1250];



%load spikes
[spikes] = bz_LoadPhy;
%get spike times OUT pulse
[status, interval] = cellfun(@(a) InIntervals(a,optoStim.timestamps),spikes.times, 'uni',false);

for iUnit = 1:length(spikes.times)
    spikes.nopulse{iUnit} = spikes.times{iUnit}(~status{iUnit});
end

spkInstPhaseStruct.times = {spikes.nopulse{iAAC}};

%define the number of frequencies with which to look at phase
%preference
nfreq   = 25; % was echt 200 of zo LK % set number of frequencies
freq    = logspace(log10(1),log10(250),nfreq); % set frequencies along a log scale
ph_bin  = linspace(-pi,pi,32); % set phases into 32 bins
ph_rate = nan(nfreq-1,32,length(spkInstPhaseStruct.times)); % reserve ph_rate space for the number of frequencies and how they are binned
%%
%hardcoded for theta
iFreqBin = 9;%1:nfreq-1 %thetabin 6.2 - 7.9

%filter LFP
filtered.data = BandpassFilter(double(d), 1250 ,[freq(iFreqBin) freq(iFreqBin+1)]); % Filter the data for the specified frequency made above
%get inst phaes and amplitude
filtered.phase  = InstPhase(filtered.data); % pull the instantaneous phase using the hilbert transform
filtered.amp    = InstAmplitude(filtered.data); % pull the instantaneous amplitude using hilbert transorm
filtered.timestamps     = ts ; % set timestamps of the filtered signal
filtered.samplingRate   = 1250; % sampling rate
filtered.filterparms.passband   = [freq(iFreqBin) freq(iFreqBin+1)]; % record passband range of the filter for that loop

spkInstPhaseStruct.times{1}     = sort(cell2mat(spkInstPhaseStruct.times')); % overbodig, want maar 1 cell - puts all cells from the recording into one variable I am not sure if this is a good idea
spkInstPhaseStruct.UID  = iAAC;

%%

%Get Power/Phase at each spike
spkInstPhaseStruct.amp = cellfun(@(X) interp1(filtered.timestamps,filtered.amp,X,'nearest'),...
    spkInstPhaseStruct.times,'uniformoutput',false);
spkInstPhaseStruct.phase = cellfun(@(X) interp1(filtered.timestamps,filtered.phase,X,'nearest'),...
    spkInstPhaseStruct.times,'uniformoutput',false);
spkInstPhaseStruct.ph_bin = ph_bin;
% histoEdges = ph_bin;
% binnedHisto = histcounts(spkInstPhaseStruct.phase{1},histoEdges);

% save([basename '.phamp.analysis.mat'],'spkInstPhaseStruct')
% histogram(binnedHisto,histoEdges)


