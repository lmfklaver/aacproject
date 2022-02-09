function [ph_portrait] = getPhasePortrait(basepath, varargin)
% This function analyzes the phase locking of each cell to the various
% frequencies wanting to be analyzed
% Rate is higher when at certain phase in certain frequencies
%
%   USAGE
%   [ph_portrait] = getPhasePortrait(basepath, <options>)
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
%   ph_portrait     - struct containing the following fields:
%   .ph_rate        - freqbin x phasebins x unit of spikerate
%   .ph_map         - struct with trans_phase and time
%   .ph_bin         - phasebins
%   .freq           -  frequencies over which phasemap is calculated
%   .nfreq          -  number of frequency bins
%
%   EXAMPLES
%   [ph_portrait] = getPhasePortrait(basepath, 'epochs', run.epochs,'saveMat',true)
%
%   TO-DO
%   - Extract few lines of code actually used from bz_PowerPhaseRatemap_sam
%   - If multishank, get ripple channel per shank in stead of one for all
%
%   HISTORY
%   Original code taken from by Sam Mckenzie
%   Edited by Kaiser Arndt (7/12/20), Lianne (Sept 2020)
%   Removed ripmod (Lianne, Dec 2020), added epochs selection (already
%   takes out of stim)


%% Parse!

if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);


p = inputParser;
addParameter(p,'basename',basename,@isstr);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'saveAs','.ph_portrait.analysis.mat');
addParameter(p,'nfreq',25,@isnumeric);
addParameter(p,'freqRange',[1 250],@isnumeric);
addParameter(p,'freqspace','log',@isstr);
addParameter(p,'muaThres',0.05,@isnumeric);
addParameter(p,'nPhaseBins',32,@isnumeric);
addParameter(p,'epochs',[],@isnumeric);


parse(p,varargin{:});
basename            = p.Results.basename;
saveMat             = p.Results.saveMat;
saveAs              = p.Results.saveAs;
nfreq               = p.Results.nfreq;
freqRange           = p.Results.freqRange;
freqspace           = p.Results.freqspace;
muaThres            = p.Results.muaThres;
nPhaseBins          = p.Results.nPhaseBins;
epochs              = p.Results.epochs;

cd(basepath)

%% Load dependencies
load([basename '.gd_eps.mat'],'gd_eps')
load([basename '.ripples.events.mat'],'ripples')
ripChan      = ripples.detectorinfo.detectionchannel;


%%%% THIS IS NECESSARY FOR PH_MAP <- SHOULD THIS BE IN HERE?
% % %% Load dependencies
% % load([basename  '.STP.mat']) % edit LK

% % %% Unpack STP
% % 
% % gd_eps              = STP.gd_eps;
% % pre_idx             = STP.pre_idx;
% % mono_con_axax       = STP.mono_con_axax;
% % post_idx            = STP.post_idx;
% % %%%%% 

%% LFP

% Get correct LFP channel (ripple channel)
xmln    = [basepath filesep basename '.xml'];
fname   = [basepath filesep basename '.lfp'];
xml     = LoadXml(xmln);

% Load in the LFP during good epochs (gd_eps) outside of stim periods
d   =[];
ts  =[];

for k = 1:size(gd_eps,1)
    dt = LoadBinary(fname,'nchannels',xml.nChannels,'channels',ripChan+1,'frequency',xml.lfpSampleRate,'start',gd_eps(k,1),'duration',diff(gd_eps(k,:)));
    d = [d;dt]; %concat all gd_eps
    ts = [ts ; [0:length(dt)-1]'/xml.lfpSampleRate + gd_eps(k,1)];
end

% Load spikes
load([basename  '.spikes.cellinfo.mat'])

%% Set Frequencies

if strcmp(freqspace,'log')
    freq = logspace(log10(freqRange(1)),log10(freqRange(2)),nfreq);
elseif strcmp(freqspace,'lin')
    freq = linspace(freqRange(1),freqRange(2),nfreq);
end

ph_bin  = linspace(-pi,pi,nPhaseBins);
ph_rate = nan(nfreq-1,nPhaseBins,length(spikes.times));

%% Loop over frequencies
for iFreq = 1:nfreq-1
    
    %filter LFP
    filter.data = BandpassFilter(double(d), 1250 ,[freq(iFreq) freq(iFreq+1)]); % Filter the data for the specified frequency made above
    
    %get inst phase and amplitude
    filtered.phase          = InstPhase(filter.data ); % inst phase using the hilbert transform
    filtered.amp            = InstAmplitude(filter.data ); % inst amp using hilbert transorm
    filtered.timestamps     = ts ;
    filtered.sampleRate     = xml.lfpSampleRate;
    filtered.filterparms.passband = [freq(iFreq) freq(iFreq+1)]; % record passband range of the filter for that loop
    spikes1.times{1}        = sort(cell2mat(spikes.times'));
    % puts all cells from the recording into one variable to calculate MUA
    % for determining powerbins to keep (Sam was not sure if this is a good
    % idea)
    
    spikes1.UID = 1;
    [PowerPhaseRatemap,~] = bz_PowerPhaseRatemap_sam(spikes1,filtered,'powernorm','none');
    
    % Changed to correct versioning with old buzcode from Sam Mckenzie,
    % because otherwise it misses a phmod field (Kaiser 5/14/20)
    
    %% Inclusion power bins (bins that pass MUA Threshold)
    % Select power bins to keep for calculating phase, bins that pass the
    % power threshold
    kp_pwr_bns = PowerPhaseRatemap.powerbins(PowerPhaseRatemap.powerskew>muaThres);
    
    % Below here needs to be fixed for looping over many recordings
    if any(kp_pwr_bns)
        kp_pwr_bns = kp_pwr_bns(1); %keeping the first value from the modulated power
        
        %convert inst amplitude to normalized amplitude (as was done in
        %bz_PowerPhaseRatemap)
        
        amp = NormToInt(log10(filtered.amp),'Z',[0 Inf],filtered.sampleRate);
        
        filtered.phase(amp<kp_pwr_bns) = nan;
        % loop over cells to get phase preference at that frequency and
        % correlation with fluctuations in that frequency band
     
        occupancy_bin = histc(filtered.phase,ph_bin); %
        
        %loop over cells
        for iUnit = 1:length(spikes.times)
           
            % first only take spikes within specified epochs
            if ~isempty(epochs)
                [status] = InIntervals(spikes.times{iUnit}, epochs);
                spikes.times{iUnit} = spikes.times{iUnit}(status);
            end
            
            %only take spikes outside of gd_eps (no stim periods)
            [status] = InIntervals(spikes.times{iUnit},gd_eps);
            if ~isempty(status)
            %smooth spike train
            k = gaussian2Dfilter([100 1],[10 1]);
            instRate = nanconvn(histc(spikes.times{iUnit}(status), filtered.timestamps(1:end)),k);
            
            %get spike phase distribution per freq. band %% nb check
            % because this looks redundant with bz_getPhasemaps
            ref  =spikes.times{iUnit}(status);
            
            spk_ph =   interp1(filtered.timestamps, filtered.phase,ref,'nearest');
           [spk_ph_cnt,spk_ph_bin] = histc(spk_ph,ph_bin);
            ph_rate(iFreq,:,iUnit) = spk_ph_cnt./occupancy_bin;
%             spk_ph_bin(iFreq,:,iUnit) = spk_ph_bin;


               %%% ASK EARL IF HE USES THE PH_MAP AT ALL %%%%
            % get spike trans (must have loaded shorttermplasticity) should
            % this be in phasemap code ?? LK
            
            %check if cell is presynaptic
% %             if ismember(iUnit,pre_idx(:,3),'rows')
% %                 
% %                 %get all post syn
% %                 
% %                 kp_post = post_idx(ismember(pre_idx(:, 3), iUnit, 'rows'), 3);
% %                 
% %                 %get post synaptic cells
% %                 for k = 1:length(kp_post)
% %                     ix = find(ismember(mono_con_axax, [iUnit kp_post(k)], 'rows'));
% %                     target = spikes.times{kp_post(k)};
% %                     
% %                     %restrict to good epochs
% %                     [status] = InIntervals(target,gd_eps);
% %                     target = target(status);
% %                     
% %                     ref = ref(spk_ph_bin>0);
% %                     spk_ph_bin = spk_ph_bin(spk_ph_bin>0);
% %                     
% %                     times = [ref;target];
% %                     groups = [ones(length(ref),1); 2*ones(length(target),1)];
% %                     
% %                     %condition over phase bin
% %                     conditions = [spk_ph_bin;ones(length(target),1)]; % why so many ones
% %                     
% %                     %run conditional CCG
% %                     [times,b] = sort(times);
% %                     groups = groups(b);
% %                     conditions = conditions(b);
% %                     duration = .2;  %hardcoded
% %                     binsize = .0008; %hardcoded
% %                     conv_wind = .015; %hardcoded
% %                     [ccg,n,t] = CCGCond(times,groups,conditions,'binsize',binsize,'duration',duration,'across_groups',false);
% %                     ccg = double(squeeze(ccg(:,1,2,:))');
% %                     
% %                     %get trans prob at every ISI
% %                     for ii = 1:size(ccg,1)
% %                         [ph_map(ix).trans_phase(iFreq,ii),prob,~] = GetTransProb(ccg(ii,:)',n(1,1,2,ii),binsize,conv_wind);
% %                         
% %                     end
% %                     
% %                     %save
% %                     ph_map(ix).time.n(iFreq,:) = squeeze(n(:,1,2,1));
% %                     ph_map(ix).time.t_lag = t;
% %                     ph_map(ix).time.bins = ph_bin;
% %                     
% %                 end
% %             end
            %%%%%%
            end
        end
    end
    
    iFreq
end
%% Condense into output

ph_portrait.ph_rate          = ph_rate;
% ph_portrait.ph_map           = ph_map;
ph_portrait.ph_bin           = ph_bin;
ph_portrait.freq             = freq;
ph_portrait.nfreq            = nfreq;

if saveMat
    save([basename saveAs],'ph_portrait')
end
%%
end
