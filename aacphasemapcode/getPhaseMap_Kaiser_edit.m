function [ph_mod] = getPhaseMap_Kaiser_edit(basepath, varargin)
% This function analyzes the phase locking of each cell to the various
% frequencies wanting to be analyzed
%
%
%   USAGE 
%   [ph_mod] = getPhaseMap_Kaiser_edit(basepath, <options>)
%
%   Dependencies: Buzcode, Englishlab\utilities, output
%   ShortTermPlasticity_axax
%   Ripple times need to be found and the channel number in the description
%   field of ".evt.rip"
%   stim times should be ".evt.ait"
%
%
%   INPUTS
%   Name-value paired inputs:
%   'basepath'      - folder in which .STP.mat and 'phmod.mat' can be found 
%                   (Required input, Default is pwd)
%   'basename'      - basefile name to load
%   'nfreq'         - How many frequencies to look at for phasemap
%                   (Default: 25)
%   'freqRange'     - (Default: [1 250])
%   'saveMat'       - 
%   'muaThres'      - MUA phase-locking threshold to verify that there is 
%                   enough power (Default: 0.05)
%   'nPhaseBins'    - (Default:32)
%   'gammaRange'
%   'thetaRange'
%   ''
%
%   OUTPUTS
%
%
%   EXAMPLES
%
%
%   NOTES
%   % ph_mod.ripmodtime is now ph_mod.ripmod.time
%   % ph_mod.ripmod is now ph_mod.ripmod.mod
%
%   TO-DO
%   - Now freqs are in log scale, make option 'lin'
%   - Extract few lines of code actually used from bz_PowerPhaseRatemap_sam
%   - Un-hardcode the frequency bins for gamma and theta
%   - If multishank, get ripple channel per shank
%
%   HISTORY
%   Code originially written by Sam Mckenzie
%   Code addapted/ edited by Kaiser Arndt (7/12/20)
%   Adapted LK 2020 sept (nb redundant code somewhere)

%% Parse! 

if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);


p = inputParser;
addParameter(p,'basename',basename,@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'nfreq',25,@isnumeric);
addParameter(p,'freqRange',[1 250],@isnumeric);
addParameter(p,'muaThres',0.05,@isnumeric);
addParameter(p,'nPhaseBins',32,@isnumeric)



parse(p,varargin{:});
basename        = p.Results.basename;
saveMat         = p.Results.saveMat;
nfreq           = p.Results.nfreq;
freqRange       = p.Results.freqRange;
muaThres        = p.Results.muaThres;
nPhaseBins      = p.Results.nPhaseBins;

cd(basepath)

%% Load dependencies
load([basename  '.STP.mat']) % edit LK

%% Unpack STP

gd_eps              = STP.gd_eps;
pre_idx             = STP.pre_idx;
mono_con_axax       = STP.mono_con_axax;
post_idx            = STP.post_idx;

%% Get correct LFP channel (ripple channel)

xmln    = [basepath filesep basename '.xml'];
fname   = [basepath filesep basename '.lfp'];
xml     = LoadXml(xmln);
ch      = ripples.detectorinfo.detectionchannel;

% Load in the LFP during good epochs (gd_eps) outside of stim periods

d =[];ts =[];
for k = 1:size(gd_eps,1)  
    dt = LoadBinary(fname,'nchannels',xml.nChannels,'channels',ch+1,'frequency',xml.lfpSampleRate,'start',gd_eps(k,1),'duration',diff(gd_eps(k,:)));  
    d = [d;dt];
    ts = [ts ; [0:length(dt)-1]'/xml.lfpSampleRate + gd_eps(k,1)];
end

% Load spikes
spikes = bz_GetSpikes('basepath',basepath);

freq = logspace(log10(freqRange(1)),log10(freqRange(2)),nfreq); % set frequencies along a log scale
ph_bin = linspace(-pi,pi,nPhaseBins); 
ph_rate = nan(nfreq-1,nPhaseBins,length(spikes.times)); 


%%
%loop over frequencies
for i = 1:nfreq-1
    
    %filter LFP
    filter.data = BandpassFilter(double(d), 1250 ,[freq(i) freq(i+1)]); % Filter the data for the specified frequency made above
    
    %get inst phase and amplitude
    filtered.phase  = InstPhase(filter.data ); % inst phase using the hilbert transform
    filtered.amp    = InstAmplitude(filter.data ); % inst amp using hilbert transorm
    filtered.timestamps     = ts ; 
    filtered.samplingRate   = xml.lfpSamplingRate;
    filtered.filterparms.passband = [freq(i) freq(i+1)]; % record passband range of the filter for that loop
    spikes1.times{1} = sort(cell2mat(spikes.times')); 
    % puts all cells from the recording into one variable to calculate MUA
    % for determining powerbins to keep (Sam was not sure if this is a good
    % idea)
    
    spikes1.UID = 1;
    [PowerPhaseRatemap,~] = bz_PowerPhaseRatemap_sam(spikes1,filtered,'powernorm','none');

    % Changed to correct versioning with old buzcode from Sam Mckenzie,
    % because otherwise it misses a phmod field (Kaiser 5/14/20)
    
    % Select power bins to keep for calculating phase, bins that pass the
    % power threshold
    kp_pwr_bns = PowerPhaseRatemap.powerbins(PowerPhaseRatemap.powerskew>muaThres);
    
    % Below here needs to be fixed for looping over many recordings
    if any(kp_pwr_bns)
        kp_pwr_bns = kp_pwr_bns(1); %keeping the first value from the modulated power
        
        %convert inst amplitude to normalized amplitude (as was done in
        %bz_PowerPhaseRatemap)
        %
        amp = NormToInt(log10(filtered.amp),'Z',[0 Inf],filtered.samplingRate);
        
        filtered.phase(amp<kp_pwr_bns) = nan;
        % loop over cells to get phase preference at that frequency and
        % correlation with fluctuations in that frequency band
        
        occ_bin = histc(filtered.phase,ph_bin);
        
        %loop over cells
        for j = 1:length(spikes.times)
            
            %only take spikes outside of gd_eps (no stim periods)
            [status] = InIntervals(spikes.times{j},gd_eps); 
            k = gaussian2Dfilter([100 1],[10 1]);
            
            %smooth spike train
            instRate = nanconvn(histc(spikes.times{j}(status), filtered.timestamps(1:end)),k);
            
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %correlate smooth spike train with inst. amplitude
%             tmp =  corrcoef(instRate,amp(1:end));
%             amp_cor(j,i) = tmp(1,2); % amp_cor is never used again in the function
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ref  =spikes.times{j}(status); 
            
            %get spike phase distribution per freq. band %% nb check
            % because this looks redundant with bz_getPhasemaps
            spk_ph =   interp1(filtered.timestamps, filtered.phase,ref,'nearest');
            [spk_ph_cnt,spk_ph_bin] = histc(spk_ph,ph_bin);
            ph_rate(i,:,j) = spk_ph_cnt./occ_bin;
            
            % get spike trans (must have loaded shorttermplasticity) should
            % this be in phasemap code ?? LK
            
            %check if cell is presynaptic
            if ismember(j,pre_idx(:,3),'rows')
                
                %get all post syn
                
                kp_post = post_idx(ismember(pre_idx(:, 3), j, 'rows'), 3);
                
                %get post synaptic cells
                for k = 1:length(kp_post)
                    ix = find(ismember(mono_con_axax, [j kp_post(k)], 'rows'));
                    target = spikes.times{kp_post(k)};
                    
                    %restrict to good epochs
                    [status] = InIntervals(target,gd_eps);
                    target = target(status);
                    
                    ref = ref(spk_ph_bin>0);
                    spk_ph_bin = spk_ph_bin(spk_ph_bin>0);
                    
                    times = [ref;target];
                    groups = [ones(length(ref),1); 2*ones(length(target),1)];
                    
                    %condition over phase bin
                    conditions = [spk_ph_bin;ones(length(target),1)]; % why so many ones
                    
                    %run conditional CCG
                    [times,b] = sort(times);
                    groups = groups(b);
                    conditions = conditions(b);
                    duration = .2;
                    binsize = .0008;
                    conv_wind = .015;
                    [ccg,n,t] = CCGCond(times,groups,conditions,'binsize',binsize,'duration',duration,'across_groups',false);
                    ccg = double(squeeze(ccg(:,1,2,:))');
                    
                    %get trans prob at every ISI
                    for ii = 1:size(ccg,1)
                        [ph_map(ix).trans_phase(i,ii),prob,~] = GetTransProb(ccg(ii,:)',n(1,1,2,ii),binsize,conv_wind);
                        
                    end
                    
                    %save
                    ph_map(ix).time.n(i,:) = squeeze(n(:,1,2,1));
                    ph_map(ix).time.t_lag = t;
                    ph_map(ix).time.bins = ph_bin;
                    
                end
            end
            
        end
        
    end
    i
    
end

%% get preferred phases for ripples and gamma

% Make this to select a bin from 'gammaRange' and 'thetaRange'
%This depends highly on how many bins! 
gammafreqbin = 18;% 39.7 - 49.95
thetafreqbin = 9; %6.3 - 7.9 Hz



ph_rate1 = ph_rate;
ph_pref_gam =[];ph_pref_theta=[];
for i = 1:size(ph_rate1,3)
    y = ph_rate1(:,:,i);
    rvect = nanmean(y(gammafreqbin,:).*exp(1i .* ph_bin)); % 
    ph_pref_gam(:,i) = atan2(imag(rvect),real(rvect));
    y = ph_rate1(:,:,i);
    rvect = nanmean(y(thetafreqbin,:).*exp(1i .* ph_bin)); %
    ph_pref_theta(:,i) = atan2(imag(rvect),real(rvect)); 
end


%% Calculate the ripple CCG and get the ripple modulation:
[ripple_ccg] = getRipCCG(basepath,varargin);
ripmod = getRipMod(basepath,'ccg',ripple_ccg);


%% Condense into output

ph_mod.ripmod           = ripmod; %with ripmod.mod and ripmod.time
ph_mod.rip_ccg          = rip_ccg;
ph_mod.rip_ccglength    = totalCCGLength;
ph_mod.ph_pref_gam      = ph_pref_gam;
ph_mod.ph_pref_theta    = ph_pref_theta;
ph_mod.ph_rate          = ph_rate;
ph_mod.ph_map           = ph_map;
ph_mod.ph_bin           = ph_bin;
ph_mod.ph_freq_theta    = freq(thetafreqbin);
ph_mod.ph_freq_gamma    = freq(gammafreqbin);
ph_mod.freq             = freq;
ph_mod.nfreq            = nfreq;



save([basename '.ph_mod.mat'], 'ph_mod')

