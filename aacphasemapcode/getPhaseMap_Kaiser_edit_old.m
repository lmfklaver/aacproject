function getPhaseMap_Kaiser_edit_old(basepath)
% This function analyzes the phase locking of each cell to the various
% frequencies wanting to be analyzed
%
%%% Dependencies %%%
% 
% Requires buzcode functions
% Requires output from functions getOptoStim and ShortTermPlasticity
% Ripple times need to be found and the channel number in the description field
% stim times should be ".evt.ait"
%
%
%%% Inputs %%%
%
% basepath (default = cd)
%
%%% To Do's %%%
%
% Make number of frequencies an input
% Make frequency range an input
% Make cell clustering function
% Make a function to run over all functions in this package (possibly leave to Lianne)
% Put plotting chunks into their own functions?
%
% Code originially written by Sam Mckenzie
% Code addapted/ edited by Kaiser Arndt (7/12/20)
%
%% Defaults
% basepath = cd;

%% Inputs

% basepath = basepath;
%% load output from short term plasticity analysis

basename = bz_BasenameFromBasepath(basepath);

load([basename  '.STP.mat']) % edit LK
% load([basename '.optmod.mat'])


%% Unpack STP

gd_eps = STP.gd_eps;
pre_idx = STP.pre_idx;
mono_con_axax_idx = STP.mono_con_axax_idx;
mono_con_axax = STP.mono_con_axax;
post_idx = STP.post_idx;
prob_uncor = STP.prob_uncor;
acg_pre = STP.acg_pre;
acg_post = STP.acg_post;
prob = STP.prob;


%% get correct LFP channel (one from ripple)

% stims isn't used in the loop possibly delete
% load stims
% fils  = getAllExtFiles(basepath,'ait',1);
% stims = LoadEvents(fils{1});


% get ripple
fils  = getAllExtFiles(basepath,'rip',1);
rip = LoadEvents(fils{1});

% pulls the channel from the ripples and loads the xml file
ch = str2num(rip.description{1}(regexp(rip.description{1},'[0-9]')));
xmln = [basepath '/' basename '.xml'];
fname = [basepath '/' basename '.lfp'];
xml = LoadXml(xmln);

% load lfp around stim
d =[];ts =[];
for k = 1:size(gd_eps,1)
    dt = LoadBinary(fname,'nchannels',xml.nChannels,'channels',ch+1,'frequency',xml.lfpSampleRate,'start',gd_eps(k,1),'duration',diff(gd_eps(k,:)));
%     lfp = bz_GetLFP(ch); % figure out why load binary gives a shorter number samples than with bz_GetLFP probably because of the gd_eps duration
%     dt = lfp.data;
    d = [d;dt];
    ts = [ts ; [0:length(dt)-1]'/1250 + gd_eps(k,1)];
end

%load spikes
spikes = bz_GetSpikes('basepath',basepath);

%define the number of frequencies with which to look at phase
%preference
nfreq = 25; % was echt 200 of zo LK % set number of frequencies
freq = logspace(log10(1),log10(250),nfreq); % set frequencies along a log scale
ph_bin = linspace(-pi,pi,32); % set phases into 32 bins
ph_rate = nan(nfreq-1,32,length(spikes.times)); % reserve ph_rate space for the number of frequencies and how they are binned
%%
%loop over frequencies
for i = 1:nfreq-1
    
    %filter LFP
    filter.data = BandpassFilter(double(d), 1250 ,[freq(i) freq(i+1)]); % Filter the data for the specified frequency made above
%     [b,a] = butter(1, [freq(i) freq(i+1)]/(1250/2), 'bandpass'); %Lianne and kaiser edit
%     filter.data = filtfilt(b,a,double(d));
    %get inst phaes and amplitude
    filtered.phase = InstPhase(filter.data ); % pull the instantaneous phase using the hilbert transform
    filtered.amp = InstAmplitude(filter.data ); % pull the instantaneous amplitude using hilbert transorm
    filtered.timestamps = ts ; % set timestamps of the filtered signal
    filtered.samplingRate  = 1250; % sampling rate
    filtered.filterparms.passband = [freq(i) freq(i+1)]; % record passband range of the filter for that loop
    spikes1.times{1} = sort(cell2mat(spikes.times')); % puts all cells from the recording into one variable I am not sure if this is a good idea
    spikes1.UID = 1;
    
    %get power threshold with MUA phase locking >.05 (arbitray
    %magnitude cut off)
    [PowerPhaseRatemap,spikebinIDs] = bz_PowerPhaseRatemap_sam(spikes1,filtered,'powernorm','none'); %... possibly optimized for Dans old recordings
    % This may have issues... just fyi...
    % Changed to correct versioning with old buzcode from Sam Mckenzie (Kaiser 5/14/20)
    
    
    kp_pwr_bns = PowerPhaseRatemap.powerbins(PowerPhaseRatemap.powerskew>.05); % setting powerbins where the frequencies are worth keeping
% below here needs to be fixed for looping over many recordings
    if any(kp_pwr_bns)
        kp_pwr_bns = kp_pwr_bns(1); %keeping the first value from the modulated
        
        %convert inst amplitude to normalized amplitude (as was done in
        %bz_PowerPhaseRatemap)
        %
        amp = NormToInt(log10(filtered.amp),'Z',[0 Inf],filtered.samplingRate);
        
        
        filtered.phase(amp<kp_pwr_bns) = nan;
        % loop over cells to get phase preference at that frequency and
        % correlation with fluctuations in that frequency band
        
        occ_bin = histc( filtered.phase,ph_bin);
        
        %loop over cells
        for j = 1:length(spikes.times)
            
            %only take spikes outside of stim periods
            [status] = InIntervals(spikes.times{j},gd_eps); %gd_eps = outside stim times
            k = gaussian2Dfilter([100 1],[10 1]);
            
            %smooth spike train
            instRate = nanconvn(histc(spikes.times{j}(status), filtered.timestamps(1:end)),k);
            
            %correlate smooth spike train with inst. amplidute
            tmp =  corrcoef(instRate,amp(1:end));
            amp_cor(j,i) = tmp(1,2); % amp_cor is never used again in the function
            ref  =spikes.times{j}(status);
            
            %get spike phase distribution per freq. band
            spk_ph =   interp1(filtered.timestamps, filtered.phase,ref,'nearest');
            [spk_ph_cnt,spk_ph_bin] = histc(spk_ph,ph_bin);
            ph_rate(i,:,j) = spk_ph_cnt./occ_bin;
            
            %get spike trans (must have loaded shorttermplasticity)
            
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

save 'phasemapcode_halfway.mat'

%% get preferred phases for ripples and gamma

% possible chunck where the issue of assigning cells to subgroups

%This depends highly on how many bins! 
gammafreqbin = 18;% 49.95
thetafreqbin = 8; %5.0049

ph_rate1 = ph_rate;
ph_pref_gam =[];ph_pref_theta=[];
for i = 1:size(ph_rate1,3)
    y = ph_rate1(:,:,i);
    rvect = nanmean(y(gammafreqbin,:).*exp(1i .* ph_bin)); % 60th frequency bin of the 200 frequencies analyzed
    ph_pref_gam(:,i) = atan2(imag(rvect),real(rvect));
    y = ph_rate1(:,:,i);
    rvect = nanmean(y(thetafreqbin,:).*exp(1i .* ph_bin)); % 35th frequency bin of the 200 frequencies analyzed
    ph_pref_theta(:,i) = atan2(imag(rvect),real(rvect)); 
end

% the indexes for theta freq and gamma freq are now hardcoded


%%



% get ripple CCGs
% reservations were all made to make this section work 
cid = [];
ripmod = [];
rip_ccg = [];
NN = [];
ix = 1;

rip = LoadEvents([basepath '/' basename '.evt.rip']);
t = rip.time(cellfun(@any,regexp(rip.description,'start')));
for j = 1:length(spikes.times)
    [status] = InIntervals(spikes.times{j},gd_eps); 
    cid = [cid;spikes.shankID(j) spikes.cluID(j)]; 
    rip_ccg(ix,:) = CrossCorr(t,spikes.times{j}(status),.005,10001)/length(t);
    NN(ix) = sum(status);
    ix = ix+1;
end


ripmod = nanmean(rip_ccg(:,5000:5025),2)./nanmean(rip_ccg(:,1:4000),2);

%% Condense into output

ph_mod.ripmod = ripmod;
ph_mod.rip_ccg = rip_ccg;
ph_mod.ph_pref_gam = ph_pref_gam;
ph_mod.ph_pref_theta = ph_pref_theta;
ph_mod.ph_rate = ph_rate;
ph_mod.ph_map = ph_map;
ph_mod.ph_bin = ph_bin;
ph_mod.freq = freq;
ph_mod.nfreq = nfreq;
ph_mod.cid = cid;

save([basename '.ph_mod.mat'], 'ph_mod')

%% Depth specifics should be in the probemap.mat files
% 
% 
% % get depth info - hard coded maps for each probe type
% 
% map1 = [...
%    5	0 ;...
% 9	20;...
% 10	40;...
% 6	60;...
% 4	80;...
% 8	100;...
% 11	120;...
% 7	140;...
% 1	0;...
% 13	20;...
% 14	40;...
% 2	60;...
% 0	80;...
% 12	100;...
% 15	120;...
% 3	140;...
% 18	0;...
% 30	20;...
% 29	40;...
% 17	60;...
% 19	80;...
% 31	100;...
% 28	120;...
% 16	140;...
% 22	0;...
% 26	20;...
% 25	40;...
% 21	60;...
% 23	80;...
% 27	100;...
% 24	120;...
% 20	140;...
% ];
% 
% 
% map2 = [...
%     16	12.65;...
% 17	37.65;...
% 18	62.65;...
% 20	87.65;...
% 21	112.65;...
% 22	137.65;...
% 31	162.65;...
% 30	187.65;...
% 29	212.65;...
% 27	237.65;...
% 26	0;...
% 25	50;...
% 24	100;...
% 28	150;...
% 23	200;...
% 19	250;...
% 12	275;...
% 8	225;...
% 3	175;...
% 7	125;...
% 6	75;...
% 5	25;...
% 4	237.65;...
% 2	212.65;...
% 1	187.65;...
% 0	162.65;...
% 9	137.65;...
% 10	112.65;...
% 11	87.65;...
% 13	62.65;...
% 14	37.65;...
% 15	12.65;...
% ];
% 
% %get depth
% celdep =[];ripleng =[];
% for i = 2:length(dirN)
%     basename = bz_BasenameFromBasepath(dirN{i});
%     spikes = bz_GetSpikes('basepath',dirN{i});
%      fils  = getAllExtFiles(dirN{i},'rip',1);
%     rip = LoadEvents(fils{1});
%     xmln = [dirN{i} filesep  basename '.xml'];
%     ch = str2num(rip.description{1}(regexp(rip.description{1},'[0-9]')));
%    
%     [xml, rxml] = LoadXml(xmln);
%     
%     %convert each shank to depth
%     
%     if length(xml.AnatGrps)==1
%          deprip =   -map2(map2(:,1)==ch,2);
%          chan = map2(:,1);
%          dep = -map2(:,2);
%          
%     else
%      deprip =   -map1(map1(:,1)==ch,2);
%          chan = map1(:,1);
%          dep = -map1(:,2);
%          
%     
%     end
%      [~,b] = ismember(spikes.maxWaveformCh,chan);
%           celdep = [celdep; dep(b)-deprip];
%     
% end
% 

